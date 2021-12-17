! (C) Copyright 2017-2020 UCAR
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
use fv3jedi_field_mod,           only: checksame, fv3jedi_field, hasfield, get_field
use fv3jedi_fields_mod,          only: fv3jedi_fields
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_kinds_mod,           only: kind_real
use wind_vt_mod,                 only: a_to_d

implicit none
private
public :: fv3jedi_state

!> Fortran derived type to hold FV3JEDI state
type, extends(fv3jedi_fields) :: fv3jedi_state
contains
  procedure, public :: add_incr
  procedure, public :: analytic_IC
  procedure, public :: set_geom_orography
end type fv3jedi_state

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine add_incr(self, geom, rhs_fields)

class(fv3jedi_state), intent(inout) :: self
type(fv3jedi_geom),   intent(inout) :: geom
type(fv3jedi_field),  intent(in)    :: rhs_fields(:)

integer :: var, i, j, k, isc, iec, jsc, jec, npz
logical :: found_neg
type(fv3jedi_field), pointer :: field_pointer

real(kind=kind_real), allocatable :: rhs_ud(:,:,:), rhs_vd(:,:,:)
real(kind=kind_real), allocatable :: rhs_delp(:,:,:)

real(kind=kind_real), pointer :: self_ua(:,:,:)
real(kind=kind_real), pointer :: self_va(:,:,:)
real(kind=kind_real), pointer :: self_ud(:,:,:)
real(kind=kind_real), pointer :: self_vd(:,:,:)
real(kind=kind_real), pointer :: self_t (:,:,:)
real(kind=kind_real), pointer :: self_pt(:,:,:)
real(kind=kind_real), pointer :: self_delp(:,:,:)
real(kind=kind_real), pointer :: self_ps(:,:,:)

real(kind=kind_real), pointer :: rhs_ua(:,:,:)
real(kind=kind_real), pointer :: rhs_va(:,:,:)
real(kind=kind_real), pointer :: rhs_t (:,:,:)
real(kind=kind_real), pointer :: rhs_pt(:,:,:)
real(kind=kind_real), pointer :: rhs_ps(:,:,:)

! Convenience
! -----------
isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec
npz = geom%npz

! Handle special cases, e.g. D grid winds to A grid winds
! -------------------------------------------------------

! Get D-Grid winds if necessary
if (hasfield(rhs_fields, 'ua')) then !A-Grid in increment
  if (.not.hasfield(rhs_fields, 'ud')) then !D-Grid not in increment
    if (self%has_field('ud')) then !D-grid in state
      allocate(rhs_ud(isc:iec  ,jsc:jec+1,1:npz))
      allocate(rhs_vd(isc:iec+1,jsc:jec  ,1:npz))
      call get_field(rhs_fields, 'ua', rhs_ua)
      call get_field(rhs_fields, 'va', rhs_va)
      call a_to_d(geom, rhs_ua, rhs_va, rhs_ud, rhs_vd) !Linear
    endif
  endif
endif


! Convert ps to delp if necessary
! -------------------------------
if (hasfield(rhs_fields, 'ps')) then !ps in increment
  if (.not.hasfield(rhs_fields, 'delp')) then !delp not in increment
    if (self%has_field('delp')) then !delp in state
      allocate(rhs_delp(isc:iec,jsc:jec,1:npz))
      call get_field(rhs_fields, 'ps', rhs_ps)
      do k = 1, npz
        rhs_delp(:,:,k) = (geom%bk(k+1)-geom%bk(k))*rhs_ps(:,:,1) !TLM
      enddo
    endif
  endif
endif

!Fields to add determined from increment
do var = 1, size(rhs_fields)

  !Winds are a special case
  if (rhs_fields(var)%fv3jedi_name == 'ua') then

    if (self%has_field('ua')) then
      call get_field(rhs_fields,  'ua', rhs_ua)
      call self%get_field('ua', self_ua)
      self_ua = self_ua + rhs_ua
    endif
    if (self%has_field('ud') .and. .not.hasfield(rhs_fields, 'ud')) then
      call self%get_field('ud', self_ud)
      self_ud = self_ud + rhs_ud
    endif

  elseif (rhs_fields(var)%fv3jedi_name == 'va') then

    if (self%has_field('va')) then
      call get_field(rhs_fields,  'va', rhs_va)
      call self%get_field('va', self_va)
      self_va = self_va + rhs_va
    endif
    if (self%has_field('vd') .and. .not.hasfield(rhs_fields, 'vd')) then
      call self%get_field('vd', self_vd)
      self_vd = self_vd + rhs_vd
    endif

  elseif (rhs_fields(var)%fv3jedi_name == 't') then

    if (self%has_field('t')) then
      call get_field(rhs_fields,  't', rhs_t)
      call self%get_field('t', self_t)
      self_t = self_t + rhs_t
    endif

    if (self%has_field('pt')) then
      call get_field(rhs_fields,  'pt', rhs_pt)
      call self%get_field('pt', self_pt)
      self_pt = self_pt + rhs_pt
    endif

  elseif (rhs_fields(var)%fv3jedi_name == 'ps') then

    if (self%has_field('ps')) then
      call get_field(rhs_fields,  'ps', rhs_ps)
      call self%get_field('ps', self_ps)
      self_ps = self_ps + rhs_ps
    endif

    if (self%has_field('delp') .and. .not. hasfield(rhs_fields, 'delp')) then
      call self%get_field('delp', self_delp)
      self_delp = self_delp + rhs_delp
    endif

  else

    !Get pointer to state
    call self%get_field(rhs_fields(var)%fv3jedi_name, field_pointer)

    !Add increment to state
    field_pointer%array = field_pointer%array + rhs_fields(var)%array

    !Nullify pointer
    nullify(field_pointer)

  endif

enddo

if (allocated(rhs_ud)) deallocate(rhs_ud)
if (allocated(rhs_vd)) deallocate(rhs_vd)
if (allocated(rhs_delp)) deallocate(rhs_delp)

!Check for negative tracers and increase to 0.0
do var = 1,self%nf

  if (self%fields(var)%tracer) then

    found_neg = .false.

    do k = 1, self%fields(var)%npz
      do j = jsc, jec
        do i = isc, iec
          if (self%fields(var)%array(i,j,k) < 0.0_kind_real) then
            found_neg = .true.
            self%fields(var)%array(i,j,k) = 0.0_kind_real
          endif
        enddo
      enddo
    enddo

    !Print message warning about negative tracer removal
    if (found_neg .and. self%f_comm%rank() == 0) print*, &
             'fv3jedi_state_mod.add_incr: Removed negative values for '&
             //trim(self%fields(var)%fv3jedi_name)

  endif

enddo

end subroutine add_incr

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
