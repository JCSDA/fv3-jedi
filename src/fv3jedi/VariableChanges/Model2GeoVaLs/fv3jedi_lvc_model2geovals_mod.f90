! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_lvc_model2geovals_mod

use iso_c_binding

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module,           only: fckit_log

use datetime_mod

use fv3jedi_constants_mod, only: constoz
use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_fieldfail_mod, only: field_fail
use fv3jedi_field_mod,     only: copy_subset, field_clen
use fv3jedi_increment_mod, only: fv3jedi_increment
use fv3jedi_kinds_mod,     only: kind_real
use fv3jedi_state_mod,     only: fv3jedi_state

use height_vt_mod
use moisture_vt_mod
use pressure_vt_mod
use surface_vt_mod
use temperature_vt_mod
use wind_vt_mod

implicit none

private
public :: fv3jedi_lvc_model2geovals

type :: fv3jedi_lvc_model2geovals
  integer :: isc, iec, jsc, jec, npz
  real(kind=kind_real), allocatable :: t (:,:,:)
  real(kind=kind_real), allocatable :: q (:,:,:)
  real(kind=kind_real), allocatable :: o3(:,:,:)
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: multiply
    procedure, public :: multiplyadjoint
end type fv3jedi_lvc_model2geovals

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, bg, fg, dummyconf)

class(fv3jedi_lvc_model2geovals), intent(inout) :: self
type(fv3jedi_geom),               intent(in)    :: geom
type(fv3jedi_state),              intent(in)    :: bg
type(fv3jedi_state),              intent(in)    :: fg
type(fckit_configuration),        intent(in)    :: dummyconf

!!! DO NOT USE CONF !!!

! Grid convenience
self%isc = geom%isc
self%iec = geom%iec
self%jsc = geom%jsc
self%jec = geom%jec
self%npz = geom%npz

! Trajectory fields
if (bg%has_field('t'     )) call bg%get_field('t'     , self%t )
if (bg%has_field('sphum' )) call bg%get_field('sphum' , self%q )
if (bg%has_field('o3mr'  )) call bg%get_field('o3mr'  , self%o3)
if (bg%has_field('o3ppmv')) call bg%get_field('o3ppmv', self%o3)

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_lvc_model2geovals), intent(inout) :: self

if (allocated(self%t )) deallocate(self%t )
if (allocated(self%q )) deallocate(self%q )
if (allocated(self%o3)) deallocate(self%o3)

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine multiply(self, geom, dxm, dxg)

class(fv3jedi_lvc_model2geovals), intent(inout) :: self
type(fv3jedi_geom),               intent(inout) :: geom
type(fv3jedi_increment),          intent(in)    :: dxm
type(fv3jedi_increment),          intent(inout) :: dxg

integer :: f, ji, jj, jk, nf2do
character(len=field_clen), allocatable :: fields_to_do_(:)
character(len=field_clen), allocatable :: fields_to_do(:)
real(kind=kind_real), pointer :: field_ptr(:,:,:)

!Winds
logical :: have_winds
real(kind=kind_real), allocatable :: ua  (:,:,:)         !A-grid wind u component
real(kind=kind_real), allocatable :: va  (:,:,:)         !A-grid wind v component
real(kind=kind_real), pointer     :: ud  (:,:,:)         !D-grid wind u component
real(kind=kind_real), pointer     :: vd  (:,:,:)         !D-grid wind v component

!Virtual temperature
logical :: have_tv
real(kind=kind_real), pointer     :: q   (:,:,:)         !Specific humidity
real(kind=kind_real), pointer     :: t   (:,:,:)         !Temperature
real(kind=kind_real), allocatable :: tv  (:,:,:)         !Virtual temperature

!Humidity mixing ratio
logical :: have_qmr
real(kind=kind_real), allocatable :: qmr (:,:,:)         !Humidity mixing ratio

!Ozone mixing ratio
logical :: have_o3
real(kind=kind_real), allocatable :: o3  (:,:,:)         !Ozone mixing ratio

!Surface pressure
logical :: have_ps
real(kind=kind_real), allocatable :: ps  (:,:,:)         !Surface pressure
real(kind=kind_real), pointer     :: delp(:,:,:)         !Pressure thickness


! Identity part of the change of fields
! -------------------------------------
call copy_subset(dxm%fields, dxg%fields, fields_to_do_)


! Winds are always special case
! -----------------------------
nf2do = 0
if (allocated(fields_to_do_)) nf2do = size(fields_to_do_)

if (dxg%has_field('ua')) then
  allocate(fields_to_do(nf2do+2))
  if (allocated(fields_to_do_)) fields_to_do(1:nf2do) = fields_to_do_
  fields_to_do(nf2do+1) = 'ua'
  fields_to_do(nf2do+2) = 'va'
else
  if (allocated(fields_to_do_)) then
    allocate(fields_to_do(nf2do))
    fields_to_do(1:nf2do) = fields_to_do_
  endif
endif


! If variable change is the identity early exit
! ---------------------------------------------
if (.not.allocated(fields_to_do)) return


! Assertion on D-Grid winds
! -------------------------
if (dxg%has_field('ud')) call abor1_ftn("GeoVaLs state should not have D-Grid winds")


! Winds
! -----
have_winds = .false.
if (dxm%has_field('ud')) then
  call dxm%get_field('ud', ud)
  call dxm%get_field('vd', vd)
  allocate(ua(self%isc:self%iec,self%jsc:self%jec,self%npz))
  allocate(va(self%isc:self%iec,self%jsc:self%jec,self%npz))
  call d2a(geom, ud, vd, ua, va)
  have_winds = .true.
elseif (dxm%has_field('ua')) then
    call dxm%get_field('ua', ua)
    call dxm%get_field('va', va)
    have_winds = .true.
endif


! Virtual temperature
! -------------------
have_tv = .false.
if (allocated(self%t) .and. allocated(self%t) .and. &
    dxm%has_field('t') .and. dxm%has_field('sphum')) then
  call dxm%get_field('t', t)
  call dxm%get_field('sphum', q)
  allocate(tv(self%isc:self%iec,self%jsc:self%jec,self%npz))
  call T_to_Tv_tl(geom, self%t, t, self%q, q, tv )
  have_tv = .true.
endif


! Humidity mixing ratio
! ---------------------
have_qmr = .false.
if (allocated(self%q) .and. dxm%has_field('sphum')) then
  call dxm%get_field('sphum', q)
  allocate(qmr(self%isc:self%iec,self%jsc:self%jec,self%npz))
  call crtm_mixratio_tl(geom, self%q, q, qmr)
  have_qmr = .true.
endif


! Ozone
! -----
have_o3 = .false.
if (allocated(self%o3) .and. dxm%has_field( 'o3mr') ) then
  call dxm%get_field('o3mr', o3)
  have_o3 = .true.
  do jk = 1, self%npz
    do jj = self%jsc, self%jec
      do ji = self%isc, self%iec
        if (self%o3(ji,jj,jk) >= 0.0_kind_real) then
          o3(ji,jj,jk) = o3(ji,jj,jk)  * constoz
        else
          o3(ji,jj,jk) = 0.0_kind_real
        endif
      enddo
    enddo
  enddo
endif
if (allocated(self%o3) .and. dxm%has_field( 'o3ppmv') ) then
  call dxm%get_field('o3ppmv', o3)
  have_o3 = .true.
  do jk = 1, self%npz
    do jj = self%jsc, self%jec
      do ji = self%isc, self%iec
        if (self%o3(ji,jj,jk) >= 0.0_kind_real) then
          continue
        else
          o3(ji,jj,jk) = 0.0_kind_real
        endif
      enddo
    enddo
  enddo
endif




! Surface pressure
! ----------------
have_ps = .false.
if (dxm%has_field( 'ps')) then
  call dxm%get_field('ps', ps)
  have_ps = .true.
elseif (dxm%has_field('delp')) then
  call dxm%get_field('delp', delp)
  allocate(ps(self%isc:self%iec,self%jsc:self%jec,1))
  ps(:,:,1) = sum(delp,3)
  have_ps = .true.
endif


! Loop over the fields not found in the input state and work through cases
! ------------------------------------------------------------------------
do f = 1, size(fields_to_do)

  call dxg%get_field(trim(fields_to_do(f)),  field_ptr)

  select case (trim(fields_to_do(f)))

  case ("ua")

    if (.not. have_winds) call field_fail(fields_to_do(f))
    field_ptr = ua

  case ("va")

    if (.not. have_winds) call field_fail(fields_to_do(f))
    field_ptr = va

  case ("tv")

    if (.not. have_tv) call field_fail(fields_to_do(f))
    field_ptr = tv

  case ("ps")

    if (.not. have_ps) call field_fail(fields_to_do(f))
    field_ptr = ps

  case ("humidity_mixing_ratio")

    if (.not. have_qmr) call field_fail(fields_to_do(f))
    field_ptr = qmr

  case ("mole_fraction_of_ozone_in_air")
    if (.not. have_o3) call field_fail(fields_to_do(f))
    field_ptr = o3

  ! Simulated but not assimilated
  case ("mass_content_of_cloud_liquid_water_in_atmosphere_layer")
  case ("mass_content_of_cloud_ice_in_atmosphere_layer")
  case ("pe")
  case ("p")

  case default

    call abor1_ftn("fv3jedi_lvc_model2geovals_mod.multiply unknown field: "//trim(fields_to_do(f)) &
                   //". Not in input field and no transform case specified.")

  end select

enddo

! Copy calendar infomation
! ------------------------
dxg%calendar_type = dxm%calendar_type
dxg%date_init = dxm%date_init

end subroutine multiply

! --------------------------------------------------------------------------------------------------

subroutine multiplyadjoint(self, geom, dxg, dxm)

class(fv3jedi_lvc_model2geovals), intent(inout) :: self
type(fv3jedi_geom),               intent(inout) :: geom
type(fv3jedi_increment),          intent(in)    :: dxg
type(fv3jedi_increment),          intent(inout) :: dxm

integer :: fg, fm, ji, jj, jk, dxg_index, num_not_copied
real(kind=kind_real), pointer :: field_ptr(:,:,:)
character(len=field_clen), allocatable :: fields_to_do(:)
character(len=field_clen) :: not_copied_(size(dxm%fields))
logical, allocatable :: field_passed(:)
integer :: noassim_index

!Winds
logical :: have_awinds, have_dwinds
integer :: ua_index, va_index
real(kind=kind_real), pointer     :: ua   (:,:,:)         !A-grid wind u component
real(kind=kind_real), pointer     :: va   (:,:,:)         !A-grid wind v component
real(kind=kind_real), allocatable :: ud   (:,:,:)         !D-grid wind u component
real(kind=kind_real), allocatable :: vd   (:,:,:)         !D-grid wind v component

!Virtual temperature
logical :: have_tv
integer :: tv_index
real(kind=kind_real), pointer     :: tv   (:,:,:)         !Virtual temperature
real(kind=kind_real), allocatable :: q_tv (:,:,:)         !Specific humidity
real(kind=kind_real), allocatable :: t_tv (:,:,:)         !Temperature
real(kind=kind_real), pointer     :: tptr (:,:,:)         !Temperature

!Humidity mixing ratio
logical :: have_qmr
integer :: qmr_index
real(kind=kind_real), pointer     :: qmr  (:,:,:)         !Virtual temperature
real(kind=kind_real), allocatable :: q_qmr(:,:,:)         !Specific humidity
real(kind=kind_real), pointer     :: qptr (:,:,:)         !Specific humidity

!Ozone mixing ratio
logical :: have_o3, have_o3_mass, have_o3_ppmv, have_mfo3
integer :: mfo3_index
real(kind=kind_real), pointer     :: mfo3 (:,:,:)         !Ozone mole frac
real(kind=kind_real), allocatable :: o3   (:,:,:)         !Ozone mixing ratio

!Surface pressure
logical :: have_ps
integer :: ps_index
real(kind=kind_real), pointer     :: ps   (:,:,:)         !Surface pressure
real(kind=kind_real), allocatable :: delp (:,:,:)         !Pressure thickness


! ! Print information
! if (geom%f_comm%rank()==0) then
!   do fg = 1, size(dxg%fields)
!     print*, "Model2GeoVaLs.multiplyAD, GeoVaLs fields IN: ", trim(dxg%fields(fg)%fv3jedi_name), minval(dxg%fields(fg)%array), maxval(dxg%fields(fg)%array)
!   enddo
!   do fm = 1, size(dxm%fields)
!     print*, "Model2GeoVaLs.multiplyAD, Model fields IN:   ", trim(dxm%fields(fm)%fv3jedi_name), minval(dxm%fields(fg)%array), maxval(dxm%fields(fg)%array)
!   enddo
! endif


! Keep track of input fields passed to output
allocate(field_passed(size(dxg%fields)))
field_passed = .false.

! Loop over model fields
num_not_copied = 0
do fm = 1, size(dxm%fields)
  ! Identity if found and not winds
  if (.not.trim(dxm%fields(fm)%fv3jedi_name) == 'ua' .and. &
      .not.trim(dxm%fields(fm)%fv3jedi_name) == 'va' .and. &
      dxg%has_field( dxm%fields(fm)%fv3jedi_name, dxg_index)) then
    call dxg%get_field(dxm%fields(fm)%fv3jedi_name, field_ptr)
    dxm%fields(fm)%array = dxm%fields(fm)%array + field_ptr
    field_passed(dxg_index) = .true.
  else
    num_not_copied = num_not_copied + 1
    not_copied_(num_not_copied) = dxm%fields(fm)%fv3jedi_name
  endif
enddo

allocate(fields_to_do(num_not_copied))
fields_to_do(1:num_not_copied) = not_copied_(1:num_not_copied)


! Winds
! -----
have_awinds = .false.
have_dwinds = .false.
if (dxg%has_field( "ua", ua_index) .and. dxg%has_field( "va", va_index)) then
  call dxg%get_field('ua', ua)
  call dxg%get_field('va', va)
  if (dxm%has_field('ud')) then
    allocate(ud(self%isc:self%iec  ,self%jsc:self%jec+1,self%npz))
    allocate(vd(self%isc:self%iec+1,self%jsc:self%jec  ,self%npz))
    ud = 0.0_kind_real
    vd = 0.0_kind_real
    call d2a_ad(geom, ud, vd, ua, va)
    have_dwinds = .true.
  elseif (dxm%has_field('ua')) then
    have_awinds = .true.
  else
    call abor1_ftn("fv3jedi_lvc_model2geovals_mod.multiplyadjoint: Winds found in GeoVaLs but"// &
                   " not in the model.")
  endif
endif


! Virtual temperature
! -------------------
have_tv = .false.
if (allocated(self%t) .and. allocated(self%t) .and. dxg%has_field('tv', tv_index)) then
  call dxg%get_field('tv', tv)
  allocate(t_tv(self%isc:self%iec,self%jsc:self%jec,self%npz))
  allocate(q_tv(self%isc:self%iec,self%jsc:self%jec,self%npz))
  t_tv = 0.0_kind_real
  q_tv = 0.0_kind_real
  call T_to_Tv_ad(geom, self%t, t_tv, self%q, q_tv, tv )
  have_tv = .true.
endif


! Humidity mixing ratio
! ---------------------
have_qmr = .false.
if (allocated(self%q) .and. dxg%has_field('humidity_mixing_ratio', qmr_index)) then
  call dxg%get_field('humidity_mixing_ratio', qmr)
  allocate(q_qmr(self%isc:self%iec,self%jsc:self%jec,self%npz))
  q_qmr = 0.0_kind_real
  call crtm_mixratio_ad(geom, self%q, q_qmr, qmr)
  have_qmr = .true.
endif


! Pressure
! --------
have_ps = .false.
if (dxg%has_field( "ps", ps_index)) then
  call dxg%get_field('ps', ps)
  allocate(delp(self%isc:self%iec,self%jsc:self%jec,self%npz))
  delp = 0.0_kind_real
  do jk = 1, self%npz
    delp(:,:,jk) = delp(:,:,jk) + ps(:,:,1)
  enddo
  have_ps = .true.
endif


! Ozone
! -----
!!!
! Ozone
! -----
have_o3 = .false.
have_o3_mass = .false.
have_o3_ppmv = .false.
have_mfo3 = .false.
have_o3_mass = dxm%has_field( 'o3mr')
have_o3_ppmv = dxm%has_field( 'o3ppmv')
have_mfo3 = dxg%has_field( "mole_fraction_of_ozone_in_air", mfo3_index)

if (allocated(self%o3) .and. have_mfo3 .and. (have_o3_mass .or. have_o3_ppmv) ) then
  call dxg%get_field('mole_fraction_of_ozone_in_air', mfo3)
  allocate(o3(self%isc:self%iec,self%jsc:self%jec,self%npz))
  o3 = mfo3
  do jk = 1, self%npz
    do jj = self%jsc, self%jec
      do ji = self%isc, self%iec
        if (self%o3(ji,jj,jk) >= 0.0_kind_real .and. have_o3_mass) then
          o3(ji,jj,jk) = o3(ji,jj,jk)  * constoz ! I think this needs to be division if you want to go from mole fraction expressed in ppmv (what is used in UFOs) to kg/kg (Model) here
        else if (have_o3_ppmv .and. self%o3(ji,jj,jk) >= 0.0_kind_real ) then
          continue
        else
          o3(ji,jj,jk) = 0.0_kind_real
        endif
      enddo
    enddo
  enddo
  have_o3 = .true.
endif

! Simulated but not assimilated
if (dxg%has_field( "mass_content_of_cloud_liquid_water_in_atmosphere_layer", noassim_index)) &
  field_passed(noassim_index) = .true.
if (dxg%has_field( "mass_content_of_cloud_ice_in_atmosphere_layer", noassim_index)) &
  field_passed(noassim_index) = .true.
if (dxg%has_field( "p", noassim_index)) &
  field_passed(noassim_index) = .true.
if (dxg%has_field( "pe", noassim_index)) &
  field_passed(noassim_index) = .true.


! ! Print information
! if (geom%f_comm%rank()==0) then
!   do fg = 1, size(dxg%fields)
!     print*, "Model2GeoVaLs.multiplyAD, GeoVaLs fields: ", trim(dxg%fields(fg)%fv3jedi_name), &
!             ", passed: ", field_passed(fg)
!   enddo
! endif


! Loop over the fields not obtainable from the input
! --------------------------------------------------
do fm = 1, size(fields_to_do)

  call dxm%get_field(trim(fields_to_do(fm)), field_ptr)

  select case(trim(fields_to_do(fm)))

  case ("ud")

    if (have_dwinds) then
      field_passed(ua_index) = .true.
      field_ptr = field_ptr + ud
    endif

  case ("vd")

    if (have_dwinds) then
      field_passed(va_index) = .true.
      field_ptr = field_ptr + vd
    endif

  case ("ua")

    if (have_awinds) then
      field_passed(ua_index) = .true.
      field_ptr = field_ptr + ua
    endif

  case ("va")

    if (have_awinds) then
      field_passed(va_index) = .true.
      field_ptr = field_ptr + va
    endif

  case ("delp")

    if (have_ps) then
      field_passed(ps_index) = .true.
      field_ptr = field_ptr + delp
    endif

  case ("o3mr")

    if (have_o3) then
      field_passed(mfo3_index) = .true.
      field_ptr = field_ptr + o3
    endif
  case ("o3ppmv")

    if (have_o3) then
      field_passed(mfo3_index) = .true.
      field_ptr = field_ptr + o3
    endif


  end select

enddo

! ! Print information
! if (geom%f_comm%rank()==0) then
!   do fg = 1, size(dxg%fields)
!     print*, "Model2GeoVaLs.multiplyAD, GeoVaLs fields: ", trim(dxg%fields(fg)%fv3jedi_name), &
!             "Passed: ", field_passed(fg)
!   enddo
! endif


! Loop over fields in input not yet assigned an output field
do fg = 1, size(dxg%fields)

  if (.not. field_passed(fg)) then

    select case(trim(dxg%fields(fg)%fv3jedi_name))

    case ("tv")

      if (.not. have_tv) call field_fail(trim(dxg%fields(fg)%fv3jedi_name))
      field_passed(tv_index) = .true.
      call dxm%get_field("t", tptr)
      call dxm%get_field("sphum", qptr)
      tptr = tptr + t_tv
      qptr = qptr + q_tv

    case ("humidity_mixing_ratio")

      if (.not. have_qmr) call field_fail(trim(dxg%fields(fg)%fv3jedi_name))
      field_passed(qmr_index) = .true.
      call dxm%get_field("sphum", qptr)
      qptr = qptr + q_qmr

    case default

      call abor1_ftn("GeoVaLs field "//trim(dxg%fields(fg)%fv3jedi_name)//" has no known link "// &
                      "to fields in model state")

    end select

  endif

enddo


! Check all fields have been linked to an output field
do fg = 1, size(dxg%fields)
  if (.not. field_passed(fg)) then
    call abor1_ftn("fv3jedi_lvc_model2geovals_mod.multiplyadjoint failed to send all geoval "// &
                   "fields to a model field")
  endif
enddo


! Copy calendar infomation
! ------------------------
dxm%calendar_type = dxg%calendar_type
dxm%date_init = dxg%date_init

! ! Print information
! if (geom%f_comm%rank()==0) then
!   do fg = 1, size(dxg%fields)
!     print*, "Model2GeoVaLs.multiplyAD, GeoVaLs fields OUT: ", trim(dxg%fields(fg)%fv3jedi_name), minval(dxg%fields(fg)%array), maxval(dxg%fields(fg)%array)
!   enddo
!   do fm = 1, size(dxm%fields)
!     print*, "Model2GeoVaLs.multiplyAD, Model fields OUT:   ", trim(dxm%fields(fm)%fv3jedi_name), minval(dxm%fields(fg)%array), maxval(dxm%fields(fg)%array)
!   enddo
! endif

end subroutine multiplyadjoint

! --------------------------------------------------------------------------------------------------

end module fv3jedi_lvc_model2geovals_mod
