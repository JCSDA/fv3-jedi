! (C) Copyright 2017-2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_state_mod

! iso
use iso_c_binding

! oops uses
use datetime_mod
use dcmip_initial_conditions_test_1_2_3, only : test1_advection_deformation, &
                                                test1_advection_hadley, test3_gravity_wave
use dcmip_initial_conditions_test_4,     only : test4_baroclinic_wave

! fckit uses
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only : log

! fv3jedi uses
use fv3jedi_field_mod,           only: fv3jedi_field
use fv3jedi_fields_mod,          only: fv3jedi_fields
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_kinds_mod,           only: kind_real

implicit none
private
public :: fv3jedi_state

!> Fortran derived type to hold FV3JEDI state
type, extends(fv3jedi_fields) :: fv3jedi_state
contains
  procedure, public :: add_increment
  procedure, public :: analytic_IC
  procedure, public :: set_geom_orography
end type fv3jedi_state

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine add_increment(self, increment)

! Arguments
class(fv3jedi_state), intent(inout) :: self
type(fv3jedi_field),  intent(in)    :: increment(:)

! Locals
integer :: f, i, j, k
logical :: found_neg
type(fv3jedi_field), pointer :: state

! Loop over the increment fields and add them to the state
do f = 1, size(increment)

  !Get pointer to state
  call self%get_field(increment(f)%short_name, state)

  !Add increment to state
  state%array = state%array + increment(f)%array

  ! Disallow tracers to become negative
  if (state%tracer) then

    found_neg = .false.

    do k = 1, state%npz
      do j = state%jsc, state%jec
        do i = state%isc, state%iec
          if (state%array(i,j,k) < 0.0_kind_real) then
            found_neg = .true.
            state%array(i,j,k) = 0.0_kind_real
          endif
        enddo
      enddo
    enddo

    !Print message warning about negative tracer removal
    if (found_neg .and. self%f_comm%rank() == 0) print*, &
      'fv3jedi_state_mod.add_incr: Removed negative values for '//trim(state%long_name)

  endif

  !Nullify pointer
  nullify(state)

enddo

end subroutine add_increment

! --------------------------------------------------------------------------------------------------
!> Analytic Initialization for the FV3 Model
!!
!! \details **analytic_IC()** initializes the FV3JEDI state and State objects using one of
!! several alternative idealized analytic models.  This is intended to facilitate testing by
!! eliminating the need to read in the initial state from a file and by providing exact expressions
!! to test interpolations.  This function is activated by setting the "analytic_init" state in the
!! "initial" or "statefile" section of the configuration file.
!!
!! Initialization options that begin with "dcmip" refer to tests defined by the multi-institutional
!! 2012 [Dynamical Core Intercomparison Project](https://earthsystealcmcog.org/projects/dcmip-2012)
!! and the associated Summer School, sponsored by NOAA, NSF, DOE, NCAR, and the University of
!! Michigan.
!!
!! Currently implemented options for analytic_init include:
!! * dcmip-test-1-1: 3D deformational flow
!! * dcmip-test-1-2: 3D Hadley-like meridional circulation
!! * dcmip-test-3-1: Non-hydrostatic gravity wave
!! * dcmip-test-4-0: Baroclinic instability
!!
!! \author M. Miesch (adapted from a pre-existing call to invent_state)
!! \date March, 2018: Created
!! \date May, 2018: Added dcmip-test-3-1
!! \date June, 2018: Added dcmip-test-4-0
!!
!! \warning This routine does not initialize the fv3jedi_interp member of the fv3jedi_state object
!!
!! \warning Though an input state file is not required for these analytic initialization routines,
!! some grid information (in particular the hybrid vertical grid coefficients ak and bk)
!! is still read in from an input file when creating the geometry object that is a required
!! member of fv3jedi_state; see c_fv3jedi_geo_setup() in fv3jedi_geom_mod.F90.
!!
!!
subroutine analytic_IC(self, geom, conf)

class(fv3jedi_state),      intent(inout) :: self    !< State
type(fv3jedi_geom),        intent(inout) :: geom    !< Geometry
type(fckit_configuration), intent(in)    :: conf    !< Configuration

integer :: i,j,k
real(kind=kind_real) :: rlat, rlon
real(kind=kind_real) :: pk,pe1,pe2,ps
real(kind=kind_real) :: u0,v0,w0,t0,phis0,ps0,rho0,hum0,q1,q2,q3,q4
character(len=:), allocatable :: method
real(kind=kind_real), pointer :: ua  (:,:,:)
real(kind=kind_real), pointer :: va  (:,:,:)
real(kind=kind_real), pointer :: t   (:,:,:)
real(kind=kind_real), pointer :: delp(:,:,:)
real(kind=kind_real), pointer :: p   (:,:,:)
real(kind=kind_real), pointer :: q   (:,:,:)
real(kind=kind_real), pointer :: qi  (:,:,:)
real(kind=kind_real), pointer :: ql  (:,:,:)
real(kind=kind_real), pointer :: o3  (:,:,:)
real(kind=kind_real), pointer :: phis(:,:,:)
real(kind=kind_real), pointer :: w   (:,:,:)

! Get method for analytic intitial condition from parameters
call conf%get_or_die("method", method)

! Messages
if (geom%f_comm%rank() == 0) then
  call log%warning("fv3jedi_state: analytic initital condition with method: "//method)
endif

! Pointers to fields
call self%get_field('ua'     , ua  )
call self%get_field('va'     , va  )
call self%get_field('t'      , t   )
call self%get_field('delp'   , delp)
call self%get_field('p'      , p   )
call self%get_field('sphum'  , q   )
call self%get_field('ice_wat', qi  )
call self%get_field('liq_wat', ql  )
call self%get_field('phis'   , phis)
call self%get_field('o3mr'   , o3  )
call self%get_field('w'      , w   )

! Initialize fields
ua   = 0.0_kind_real
va   = 0.0_kind_real
t    = 0.0_kind_real
delp = 0.0_kind_real
p    = 0.0_kind_real
q    = 0.0_kind_real
qi   = 0.0_kind_real
ql   = 0.0_kind_real
phis = 0.0_kind_real
o3   = 0.0_kind_real
w    = 0.0_kind_real

select case (method)

  case ("dcmip-test-1-1")

    do i = geom%isc,geom%iec
      do j = geom%jsc,geom%jec

        rlat = geom%grid_lat(i,j)
        rlon = geom%grid_lon(i,j)

        ! Call the routine first just to get the surface pressure
        call test1_advection_deformation(rlon,rlat,pk,0.d0,1,u0,v0,w0,t0,phis0,ps,rho0,hum0,q1,q2,&
                                         q3,q4)

        phis(i,j,1) = phis0

        ! Now loop over all levels
        do k = 1, geom%npz

          pe1 = geom%ak(k) + geom%bk(k)*ps
          pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
          pk = 0.5_kind_real * (pe1+pe2)
          call test1_advection_deformation(rlon,rlat,pk,0.d0,0,u0,v0,w0,t0,phis0,ps0,rho0,hum0,q1,&
                                           q2,q3,q4)

          ua(i,j,k)   = u0
          va(i,j,k)   = v0
          t(i,j,k)    = t0
          delp(i,j,k) = pe2-pe1
          p(i,j,k)    = pk
          q (i,j,k)   = hum0
          qi(i,j,k)   = q1
          ql(i,j,k)   = q2
          o3(i,j,k)   = q3
          w(i,j,k)    = w0

        enddo
      enddo
    enddo

  case ("dcmip-test-1-2")

    do i = geom%isc,geom%iec
      do j = geom%jsc,geom%jec

        rlat = geom%grid_lat(i,j)
        rlon = geom%grid_lon(i,j)

        ! Call the routine first just to get the surface pressure
        call test1_advection_hadley(rlon,rlat,pk,0.d0,1,u0,v0,w0,t0,phis0,ps,rho0,hum0,q1)

        phis(i,j,1) = phis0

        ! Now loop over all levels
        do k = 1, geom%npz

          pe1 = geom%ak(k) + geom%bk(k)*ps
          pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
          pk = 0.5_kind_real * (pe1+pe2)
          call test1_advection_hadley(rlon,rlat,pk,0.d0,0,u0,v0,w0,t0,phis0,ps,rho0,hum0,q1)

          ua(i,j,k)   = u0
          va(i,j,k)   = v0
          t(i,j,k)    = t0
          delp(i,j,k) = pe2-pe1
          p(i,j,k)    = pk
          q(i,j,k)    = hum0
          qi(i,j,k)   = q1
          w(i,j,k)    = w0

        enddo
      enddo
    enddo

  case ("dcmip-test-3-1")

    do i = geom%isc,geom%iec
      do j = geom%jsc,geom%jec

        rlat = geom%grid_lat(i,j)
        rlon = geom%grid_lon(i,j)

        ! Call the routine first just to get the surface pressure
        call test3_gravity_wave(rlon,rlat,pk,0.d0,1,u0,v0,w0,t0,phis0,ps,rho0,hum0)

        phis(i,j,1) = phis0

        ! Now loop over all levels
        do k = 1, geom%npz

          pe1 = geom%ak(k) + geom%bk(k)*ps
          pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
          pk = 0.5_kind_real * (pe1+pe2)
          call test3_gravity_wave(rlon,rlat,pk,0.d0,0,u0,v0,w0,t0,phis0,ps,rho0,hum0)

          ua(i,j,k)   = u0
          va(i,j,k)   = v0
          t(i,j,k)    = t0
          delp(i,j,k) = pe2-pe1
          p(i,j,k)    = pk
          q(i,j,k)    = hum0
          w(i,j,k)    = w0

        enddo
      enddo
    enddo

  case ("dcmip-test-4-0")

    do i = geom%isc,geom%iec
      do j = geom%jsc,geom%jec

        rlat = geom%grid_lat(i,j)
        rlon = geom%grid_lon(i,j)

        ! Call the routine first just to get the surface pressure
        call test4_baroclinic_wave(0,1.0_kind_real,rlon,rlat,pk,0.d0,1,u0,v0,w0,t0,phis0,ps,rho0,&
                                   hum0,q1,q2)

        phis(i,j,1) = phis0

        ! Now loop over all levels
        do k = 1, geom%npz

          pe1 = geom%ak(k) + geom%bk(k)*ps
          pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
          pk = 0.5_kind_real * (pe1+pe2)
          call test4_baroclinic_wave(0,1.0_kind_real,rlon,rlat,pk,0.d0,0,u0,v0,w0,t0,phis0,ps,rho0,&
                                     hum0,q1,q2)

          ua(i,j,k)   = u0
          va(i,j,k)   = v0
          t(i,j,k)    = t0
          delp(i,j,k) = pe2-pe1
          p(i,j,k)    = pk
          q(i,j,k)    = hum0
          w(i,j,k)    = w0

        enddo
      enddo
    enddo

  case default

    call abor1_ftn("fv3jedi_state analytic_ic: provide analytic initial condition method")

end select

end subroutine analytic_ic

! --------------------------------------------------------------------------------------------------

subroutine set_geom_orography(self, geom)

  class(fv3jedi_state),      intent(in) :: self    !< State
  type(fv3jedi_geom),        intent(inout) :: geom    !< Geometry

  real(kind=kind_real), pointer :: orog(:,:,:)

  call self%get_field('orog_filt', orog)
  geom%orography = orog

end subroutine set_geom_orography

! --------------------------------------------------------------------------------------------------

end module fv3jedi_state_mod
