! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_state_mod

! iso
use iso_c_binding

! oops uses
use datetime_mod
use oops_variables_mod

! fckit uses
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only : log

! fv3jedi uses
use fields_metadata_mod,         only: field_metadata
use fv3jedi_field_mod
use fv3jedi_constants_mod,       only: rad2deg
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_increment_utils_mod, only: fv3jedi_increment
use fv3jedi_interpolation_mod,   only: field2field_interp
use fv3jedi_kinds_mod,           only: kind_real
use fv3jedi_io_gfs_mod,          only: fv3jedi_io_gfs
use fv3jedi_io_geos_mod,         only: fv3jedi_io_geos
use fv3jedi_io_latlon_mod,       only: fv3jedi_llgeom
use fv3jedi_state_utils_mod,     only: fv3jedi_state

use wind_vt_mod, only: a2d

implicit none
private
public :: fv3jedi_state, create, delete, zeros, copy, axpy, add_incr, read_file, write_file, rms, &
          change_resol, analytic_IC, getminmaxrms, &
          fv3jedi_state_serialize, fv3jedi_state_deserialize

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, vars)

implicit none
type(fv3jedi_state),  intent(inout) :: self
type(fv3jedi_geom),   intent(in)    :: geom
type(oops_variables), intent(in)    :: vars

integer :: var, fc
type(field_metadata) :: fmd

! Allocate fields structure
! -------------------------
self%nf = vars%nvars()
allocate(self%fields(self%nf))

! Loop through and allocate main state fields
! -------------------------------------------
fc = 0
do var = 1, vars%nvars()

  fmd = geom%fields%get_field(trim(vars%variable(var)))

  fc=fc+1;
  call self%fields(fc)%allocate_field(geom%isc, geom%iec, geom%jsc, geom%jec, &
                                      fmd%levels, &
                                      short_name   = trim(fmd%field_io_name), &
                                      long_name    = trim(fmd%long_name), &
                                      fv3jedi_name = trim(fmd%field_name), &
                                      units        = fmd%units, &
                                      io_file      = trim(fmd%io_file), &
                                      space        = trim(fmd%space), &
                                      staggerloc   = trim(fmd%stagger_loc), &
                                      tracer       = fmd%tracer, &
                                      integerfield = trim(fmd%array_kind)=='integer')

enddo

if (fc .ne. self%nf) &
call abor1_ftn("fv3jedi_state_mod create: fc does not equal self%nf")

self%hydrostatic = .true.
if (has_field(self%fields, 'delz') .and. has_field(self%fields, 'w')) self%hydrostatic = .false.

! Initialize all arrays to zero
call zeros(self)

! Copy some geometry for convenience
self%isc    = geom%isc
self%iec    = geom%iec
self%jsc    = geom%jsc
self%jec    = geom%jec
self%npx    = geom%npx
self%npy    = geom%npy
self%npz    = geom%npz
self%ntile  = geom%ntile
self%ntiles = geom%ntiles

! Pointer to fv3jedi communicator
self%f_comm = geom%f_comm

! Check winds
if (has_field(self%fields, 'ua') .and. .not.has_field(self%fields, 'va')) &
call abor1_ftn("fv3jedi_state_mod create: found A-Grid u but not v")
if (.not.has_field(self%fields, 'ua') .and. has_field(self%fields, 'va')) &
call abor1_ftn("fv3jedi_state_mod create: found A-Grid v but not u")
if (has_field(self%fields, 'ud') .and. .not.has_field(self%fields, 'vd')) &
call abor1_ftn("fv3jedi_state_mod create: found D-Grid u but not v")
if (.not.has_field(self%fields, 'ud') .and. has_field(self%fields, 'vd')) &
call abor1_ftn("fv3jedi_state_mod create: found D-Grid v but not u")

!Check User's choice of ozone variables.
if (has_field(self%fields, 'o3mr') .and. has_field(self%fields, 'o3ppmv')) &
call abor1_ftn("fv3jedi_state_mod create: found both o3mr and o3ppmv there can only be o3 in kgkg-1 or ppmv")


self%have_agrid = .false.
self%have_dgrid = .false.
if (has_field(self%fields, 'ua')) self%have_agrid = .true.
if (has_field(self%fields, 'ud')) self%have_dgrid = .true.

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_state), intent(inout) :: self
integer :: var

do var = 1, self%nf
  call self%fields(var)%deallocate_field()
enddo
deallocate(self%fields)

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine zeros(self)

implicit none
type(fv3jedi_state), intent(inout) :: self
integer :: var

do var = 1, self%nf
  self%fields(var)%array = 0.0_kind_real
enddo

end subroutine zeros

! --------------------------------------------------------------------------------------------------

subroutine copy(self,rhs)

implicit none
type(fv3jedi_state), intent(inout) :: self
type(fv3jedi_state), intent(in)    :: rhs

integer :: var

call checksame(self%fields,rhs%fields,"fv3jedi_state_mod.copy")

do var = 1, self%nf
  self%fields(var) = rhs%fields(var)
enddo

self%calendar_type = rhs%calendar_type
self%date_init = rhs%date_init

end subroutine copy

! --------------------------------------------------------------------------------------------------

subroutine axpy(self,zz,rhs)

implicit none
type(fv3jedi_state),  intent(inout) :: self
real(kind=kind_real), intent(in)    :: zz
type(fv3jedi_state),  intent(in)    :: rhs

integer :: var

call checksame(self%fields,rhs%fields,"fv3jedi_state_mod.axpy")

do var = 1, self%nf
  self%fields(var)%array = self%fields(var)%array + zz * rhs%fields(var)%array
enddo

end subroutine axpy

! --------------------------------------------------------------------------------------------------

subroutine add_incr(geom,self,rhs)

implicit none
type(fv3jedi_geom),      intent(inout) :: geom
type(fv3jedi_state),     intent(inout) :: self
type(fv3jedi_increment), intent(in)    :: rhs

integer :: var,i,j,k
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

! Handle special cases, e.g. D grid winds to A grid winds
! -------------------------------------------------------

! Get D-Grid winds if necessary
if (rhs%have_agrid) then !A-Grid in increment
  if (.not.rhs%have_dgrid) then !D-Grid not in increment
    if (self%have_dgrid) then !D-grid in state
      allocate(rhs_ud(rhs%isc:rhs%iec  ,rhs%jsc:rhs%jec+1,1:rhs%npz))
      allocate(rhs_vd(rhs%isc:rhs%iec+1,rhs%jsc:rhs%jec  ,1:rhs%npz))
      call pointer_field_array(rhs%fields, 'ua', rhs_ua)
      call pointer_field_array(rhs%fields, 'va', rhs_va)
      call a2d(geom, rhs_ua, rhs_va, rhs_ud, rhs_vd) !Linear
    endif
  endif
endif


! Convert ps to delp if necessary
! -------------------------------
if (has_field(rhs%fields, 'ps')) then !ps in increment
  if (.not.has_field(rhs%fields, 'delp')) then !delp not in increment
    if (has_field(self%fields, 'delp')) then !delp in state
      allocate(rhs_delp(rhs%isc:rhs%iec,rhs%jsc:rhs%jec,1:rhs%npz))
      call pointer_field_array(rhs%fields, 'ps', rhs_ps)
      do k = 1,rhs%npz
        rhs_delp(:,:,k) = (geom%bk(k+1)-geom%bk(k))*rhs_ps(:,:,1) !TLM
      enddo
    endif
  endif
endif

!Fields to add determined from increment
do var = 1,rhs%nf

  !Winds are a special case
  if (rhs%fields(var)%fv3jedi_name == 'ua') then

    if (has_field(self%fields, 'ua')) then
      call pointer_field_array(rhs%fields,  'ua', rhs_ua)
      call pointer_field_array(self%fields, 'ua', self_ua)
      self_ua = self_ua + rhs_ua
    endif
    if (has_field(self%fields, 'ud') .and. .not.has_field(rhs%fields, 'ud')) then
      call pointer_field_array(self%fields, 'ud', self_ud)
      self_ud = self_ud + rhs_ud
    endif

  elseif (rhs%fields(var)%fv3jedi_name == 'va') then

    if (has_field(self%fields, 'va')) then
      call pointer_field_array(rhs%fields,  'va', rhs_va)
      call pointer_field_array(self%fields, 'va', self_va)
      self_va = self_va + rhs_va
    endif
    if (has_field(self%fields, 'vd') .and. .not.has_field(rhs%fields, 'vd')) then
      call pointer_field_array(self%fields, 'vd', self_vd)
      self_vd = self_vd + rhs_vd
    endif

  elseif (rhs%fields(var)%fv3jedi_name == 't') then

    if (has_field(self%fields, 't')) then
      call pointer_field_array(rhs%fields,  't', rhs_t)
      call pointer_field_array(self%fields, 't', self_t)
      self_t = self_t + rhs_t
    endif

    if (has_field(self%fields, 'pt')) then
      call pointer_field_array(rhs%fields,  'pt', rhs_pt)
      call pointer_field_array(self%fields, 'pt', self_pt)
      self_pt = self_pt + rhs_pt
    endif

  elseif (rhs%fields(var)%fv3jedi_name == 'ps') then

    if (has_field(self%fields, 'ps')) then
      call pointer_field_array(rhs%fields,  'ps', rhs_ps)
      call pointer_field_array(self%fields, 'ps', self_ps)
      self_ps = self_ps + rhs_ps
    endif

    if (has_field(self%fields, 'delp') .and. .not. has_field(rhs%fields, 'delp')) then
      call pointer_field_array(self%fields, 'delp', self_delp)
      self_delp = self_delp + rhs_delp
    endif

  else

    !Get pointer to state
    call pointer_field(self%fields, rhs%fields(var)%fv3jedi_name, field_pointer)

    !Add increment to state
    field_pointer%array = field_pointer%array + rhs%fields(var)%array

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

    do k = 1,self%fields(var)%npz
      do j = geom%jsc,geom%jec
        do i = geom%isc,geom%iec
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

subroutine change_resol(self,geom,rhs,geom_rhs)

implicit none
type(fv3jedi_state), intent(inout) :: self
type(fv3jedi_geom),  intent(inout) :: geom
type(fv3jedi_state), intent(in)    :: rhs
type(fv3jedi_geom),  intent(inout) :: geom_rhs

! Interpolation
integer :: var
type(field2field_interp) :: interp
logical :: integer_interp = .false.

call checksame(self%fields,rhs%fields,"fv3jedi_state_mod.change_resol")

if ((rhs%iec-rhs%isc+1)-(self%iec-self%isc+1) == 0) then

  call copy(self, rhs)

else

  ! Check if integer interp needed
  do var = 1, self%nf
    if (rhs%fields(var)%integerfield) integer_interp = .true.
  enddo

  call interp%create(geom%interp_method, integer_interp, geom_rhs, geom)
  call interp%apply(self%nf, geom_rhs, rhs%fields, geom, self%fields)
  call interp%delete()

  self%calendar_type = rhs%calendar_type
  self%date_init = rhs%date_init

endif

end subroutine change_resol

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
!! and the associated Summer School, sponsored by NOAA, NSF, DOE, NCAR, and the University of Michigan.
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
subroutine analytic_IC(self, geom, c_conf, vdate)

  use dcmip_initial_conditions_test_1_2_3, only : test1_advection_deformation, &
                                                  test1_advection_hadley, test3_gravity_wave
  use dcmip_initial_conditions_test_4,     only : test4_baroclinic_wave

  implicit none

  type(fv3jedi_state), intent(inout) :: self    !< State
  type(fv3jedi_geom),  intent(inout) :: geom    !< Geometry
  type(c_ptr), intent(in)            :: c_conf  !< Configuration
  type(datetime), intent(inout)      :: vdate   !< DateTime

  character(len=30) :: IC
  character(len=20) :: sdate
  character(len=1024) :: buf
  Integer :: i,j,k
  real(kind=kind_real) :: rlat, rlon
  real(kind=kind_real) :: pk,pe1,pe2,ps
  real(kind=kind_real) :: u0,v0,w0,t0,phis0,ps0,rho0,hum0,q1,q2,q3,q4

  type(fckit_configuration) :: f_conf
  character(len=:), allocatable :: str

  real(kind=kind_real), pointer :: ud  (:,:,:)
  real(kind=kind_real), pointer :: vd  (:,:,:)
  real(kind=kind_real), pointer :: t   (:,:,:)
  real(kind=kind_real), pointer :: delp(:,:,:)
  real(kind=kind_real), pointer :: q   (:,:,:)
  real(kind=kind_real), pointer :: qi  (:,:,:)
  real(kind=kind_real), pointer :: ql  (:,:,:)
  real(kind=kind_real), pointer :: o3  (:,:,:)
  real(kind=kind_real), pointer :: phis(:,:,:)
  real(kind=kind_real), pointer :: w   (:,:,:)
  real(kind=kind_real), pointer :: delz(:,:,:)

  ! Fortran configuration
  ! ---------------------
  f_conf = fckit_configuration(c_conf)

  If (f_conf%has("analytic_init")) Then
     call f_conf%get_or_die("analytic_init",str)
     IC = str
     deallocate(str)
  EndIf

  call log%warning("fv3jedi_state:analytic_init: "//IC)
  call f_conf%get_or_die("date",str)
  sdate = str
  deallocate(str)
  WRITE(buf,*) 'validity date is: '//sdate
  call log%info(buf)
  call datetime_set(sdate, vdate)

  !Pointers to states
  call pointer_field_array(self%fields, 'ud'     , ud  )
  call pointer_field_array(self%fields, 'vd'     , vd  )
  call pointer_field_array(self%fields, 't'      , t   )
  call pointer_field_array(self%fields, 'delp'   , delp)
  call pointer_field_array(self%fields, 'sphum'  , q   )
  call pointer_field_array(self%fields, 'ice_wat', qi  )
  call pointer_field_array(self%fields, 'liq_wat', ql  )
  if ( has_field(self%fields, 'o3mr') ) call pointer_field_array(self%fields, 'o3mr'   , o3  )
  if ( has_field(self%fields, 'o3ppmv') ) call pointer_field_array(self%fields, 'o3ppmv'   , o3  )
  call pointer_field_array(self%fields, 'phis'   , phis)
  if (.not.self%hydrostatic) then
    call pointer_field_array(self%fields, 'w'   , w   )
    call pointer_field_array(self%fields, 'delz', delz)
  endif

  int_option: Select Case (IC)

     Case ("dcmip-test-1-1")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test1_advection_deformation(rlon,rlat,pk,0.d0,1,u0,v0,w0,t0,&
                                               phis0,ps,rho0,hum0,q1,q2,q3,q4)

              phis(i,j,1) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test1_advection_deformation(rlon,rlat,pk,0.d0,0,u0,v0,w0,t0,&
                                                  phis0,ps0,rho0,hum0,q1,q2,q3,q4)

                 ud(i,j,k) = u0 !ATTN Not going to necessary keep a-grid winds, u can be either a or d grid
                 vd(i,j,k) = v0 !so this needs to be generic. You cannot drive the model with A grid winds
                 If (.not.self%hydrostatic) w(i,j,k) = w0
                 t(i,j,k) = t0
                 delp(i,j,k) = pe2-pe1
                 q (i,j,k) = hum0
                 qi(i,j,k) = q1
                 ql(i,j,k) = q2
                 o3(i,j,k) = q3

              enddo
           enddo
        enddo

     Case ("dcmip-test-1-2")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test1_advection_hadley(rlon,rlat,pk,0.d0,1,u0,v0,w0,&
                                          t0,phis0,ps,rho0,hum0,q1)

              phis(i,j,1) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test1_advection_hadley(rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                                             t0,phis0,ps,rho0,hum0,q1)

                 ud(i,j,k) = u0 !ATTN comment above
                 vd(i,j,k) = v0
                 If (.not.self%hydrostatic) w(i,j,k) = w0
                 t(i,j,k) = t0
                 delp(i,j,k) = pe2-pe1
                 q(i,j,k) = hum0
                 qi(i,j,k) = q1

              enddo
           enddo
        enddo

     Case ("dcmip-test-3-1")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test3_gravity_wave(rlon,rlat,pk,0.d0,1,u0,v0,w0,&
                                      t0,phis0,ps,rho0,hum0)

              phis(i,j,1) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test3_gravity_wave(rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                                         t0,phis0,ps,rho0,hum0)

                 ud(i,j,k) = u0 !ATTN comment above
                 vd(i,j,k) = v0
                 If (.not.self%hydrostatic) w(i,j,k) = w0
                 t(i,j,k) = t0
                 delp(i,j,k) = pe2-pe1
                 q(i,j,k) = hum0

              enddo
           enddo
        enddo

     Case ("dcmip-test-4-0")

        do i = geom%isc,geom%iec
           do j = geom%jsc,geom%jec
              rlat = geom%grid_lat(i,j)
              rlon = geom%grid_lon(i,j)

              ! Call the routine first just to get the surface pressure
              Call test4_baroclinic_wave(0,1.0_kind_real,rlon,rlat,pk,0.d0,1,u0,v0,w0,&
                                         t0,phis0,ps,rho0,hum0,q1,q2)

              phis(i,j,1) = phis0

              ! Now loop over all levels
              do k = 1, geom%npz

                 pe1 = geom%ak(k) + geom%bk(k)*ps
                 pe2 = geom%ak(k+1) + geom%bk(k+1)*ps
                 pk = 0.5_kind_real * (pe1+pe2)
                 Call test4_baroclinic_wave(0,1.0_kind_real,rlon,rlat,pk,0.d0,0,u0,v0,w0,&
                                         t0,phis0,ps,rho0,hum0,q1,q2)

                 ud(i,j,k) = u0 !ATTN comment above
                 vd(i,j,k) = v0
                 If (.not.self%hydrostatic) w(i,j,k) = w0
                 t(i,j,k) = t0
                 delp(i,j,k) = pe2-pe1
                 q(i,j,k) = hum0

              enddo
           enddo
        enddo

     Case Default

        call abor1_ftn("fv3jedi_state analytic_IC: provide analytic_init")

     End Select int_option

end subroutine analytic_IC

! --------------------------------------------------------------------------------------------------

subroutine read_file(geom, self, c_conf, vdate)
use string_utils

implicit none

type(fv3jedi_geom),  intent(inout) :: geom     !< Geometry
type(fv3jedi_state), intent(inout) :: self     !< State
type(c_ptr),         intent(in)    :: c_conf   !< Configuration
type(datetime),      intent(inout) :: vdate    !< DateTime

type(fv3jedi_io_gfs)  :: gfs
type(fv3jedi_io_geos) :: geos

character(len=10) :: filetype
type(fckit_configuration) :: f_conf
character(len=:), allocatable :: str


! Fortran configuration
! ---------------------
f_conf = fckit_configuration(c_conf)

call f_conf%get_or_die("filetype",str)
filetype = str
deallocate(str)


if (trim(filetype) == 'gfs') then

  call gfs%setup_conf(f_conf)
  call gfs%read_meta(geom, vdate, self%calendar_type, self%date_init)
  call gfs%read_fields(geom, self%fields)

elseif (trim(filetype) == 'geos') then

  call geos%setup_conf(geom, f_conf)
  call geos%read_meta(geom, vdate, self%calendar_type, self%date_init, self%fields)
  call geos%read_fields(geom, self%fields)
  call geos%delete()

else

  call abor1_ftn("fv3jedi_state_mod.read: restart type not supported")

endif

end subroutine read_file

! --------------------------------------------------------------------------------------------------

subroutine write_file(geom, self, c_conf, vdate)

  implicit none

  type(fv3jedi_geom),  intent(inout) :: geom     !< Geometry
  type(fv3jedi_state), intent(inout) :: self     !< State
  type(c_ptr),         intent(in)    :: c_conf   !< Configuration
  type(datetime),      intent(inout) :: vdate    !< DateTime

  type(fv3jedi_io_gfs)  :: gfs
  type(fv3jedi_io_geos) :: geos
  type(fv3jedi_llgeom)  :: latlon

  character(len=10) :: filetype
  type(fckit_configuration) :: f_conf
  character(len=:), allocatable :: str


  ! Fortran configuration
  ! ---------------------
  f_conf = fckit_configuration(c_conf)


  call f_conf%get_or_die("filetype",str)
  filetype = str
  deallocate(str)

  if (trim(filetype) == 'gfs') then

    call gfs%setup_conf(f_conf)
    call gfs%setup_date(vdate)
    call gfs%write(geom, self%fields, vdate, self%calendar_type, self%date_init)

  elseif (trim(filetype) == 'geos') then

    call geos%setup_conf(geom, f_conf)
    call geos%setup_date(vdate)
    call geos%write(geom, self%fields, vdate)
    call geos%delete()

  elseif (trim(filetype) == 'latlon') then

    call latlon%setup_conf(geom)
    call latlon%setup_date(vdate)
    call latlon%write(geom, self%fields, f_conf, vdate)

  else

     call abor1_ftn("fv3jedi_state_mod.write: restart type not supported")

  endif

end subroutine write_file

! --------------------------------------------------------------------------------------------------

subroutine getminmaxrms(self, field_num, field_name, minmaxrms)

implicit none
type(fv3jedi_state),  intent(in)    :: self
integer,              intent(in)    :: field_num
character(len=*),     intent(inout) :: field_name
real(kind=kind_real), intent(inout) :: minmaxrms(3)

call field_getminmaxrms(self%fields, field_num, field_name, minmaxrms, self%f_comm)

end subroutine getminmaxrms

! --------------------------------------------------------------------------------------------------

subroutine rms(self, prms)

implicit none
type(fv3jedi_state),  intent(in)  :: self
real(kind=kind_real), intent(out) :: prms

call fields_rms(self%nf, self%fields, prms, self%f_comm)

end subroutine rms

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_state_serialize(self,vsize,vect_inc)

implicit none

! Passed variables
type(fv3jedi_state),intent(in) :: self           !< State
integer,intent(in) :: vsize                      !< Size
real(kind_real),intent(out) :: vect_inc(vsize)   !< Vector

! Local variables
integer :: ind, var, i, j, k

! Initialize
ind = 0

! Copy
do var = 1, self%nf
  do k = 1,self%fields(var)%npz
    do j = self%fields(var)%jsc,self%fields(var)%jec
      do i = self%fields(var)%isc,self%fields(var)%iec
        ind = ind + 1
        vect_inc(ind) = self%fields(var)%array(i, j, k)
      enddo
    enddo
  enddo
enddo

end subroutine fv3jedi_state_serialize

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_deserialize(self,vsize,vect_inc,index)

implicit none

! Passed variables
type(fv3jedi_state),intent(inout) :: self             !< State
integer,intent(in) :: vsize                           !< Size
real(kind_real),intent(in) :: vect_inc(vsize)         !< Vector
integer,intent(inout) :: index                        !< Index

! Local variables
integer :: ind, var, i, j, k

! Copy
do var = 1, self%nf
  do k = 1,self%fields(var)%npz
    do j = self%fields(var)%jsc,self%fields(var)%jec
      do i = self%fields(var)%isc,self%fields(var)%iec
        self%fields(var)%array(i, j, k) = vect_inc(index + 1)
        index = index + 1
      enddo
    enddo
  enddo
enddo

end subroutine fv3jedi_state_deserialize

! ------------------------------------------------------------------------------

end module fv3jedi_state_mod
