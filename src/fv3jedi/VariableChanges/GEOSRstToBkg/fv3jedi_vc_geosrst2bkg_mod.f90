! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_vc_geosrst2bkg_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use datetime_mod

use fckit_log_module, only : fckit_log

use fv3jedi_kinds_mod,   only: kind_real
use fv3jedi_geom_mod,    only: fv3jedi_geom
use fv3jedi_state_mod,   only: fv3jedi_state

use fv3jedi_field_mod, only: copy_subset, has_field, pointer_field_array

use wind_vt_mod, only: a2d, d2a
use temperature_vt_mod, only: pt_to_t, t_to_pt
use pressure_vt_mod, only: pe_to_delp, delp_to_pe, pe_to_pk, ps_to_pe
use moisture_vt_mod, only: q4_to_q2, q2_to_q4

implicit none
private

public :: fv3jedi_vc_geosrst2bkg
public :: create
public :: delete
public :: changevar
public :: changevarinverse

type :: fv3jedi_vc_geosrst2bkg
 logical :: do_wind
 logical :: do_temp
 logical :: do_pres
 logical :: do_clds
 character(len=4) :: pres_var
end type fv3jedi_vc_geosrst2bkg

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, conf)

implicit none
type(fv3jedi_vc_geosrst2bkg), intent(inout) :: self
type(fv3jedi_geom),           intent(inout) :: geom
type(fckit_configuration),    intent(in)    :: conf

character(len=:), allocatable :: str

! Select which variables to transform
! -----------------------------------

if( .not. conf%get('do_wind', self%do_wind) ) then
  self%do_wind = .true.
endif

if( .not. conf%get('do_temperature', self%do_temp) ) then
  self%do_temp = .true.
endif

if( .not. conf%get('do_pressure', self%do_pres) ) then
  self%do_pres = .true.
endif

if( .not. conf%get('do_clouds', self%do_clds) ) then
  self%do_clds = .true.
endif

if( .not. conf%has('pres_var') ) then
  self%pres_var = 'delp' ! Default to ps
else
  call conf%get_or_die("pres_var",str)
  self%pres_var = str
  deallocate(str)
endif

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_vc_geosrst2bkg), intent(inout) :: self

end subroutine delete

! ------------------------------------------------------------------------------

subroutine changevar(self,geom,xr,xb)

implicit none
type(fv3jedi_vc_geosrst2bkg), intent(inout) :: self
type(fv3jedi_geom),       intent(inout) :: geom
type(fv3jedi_state),      intent(in)    :: xr
type(fv3jedi_state),      intent(inout) :: xb

logical :: have_fractions

! Poitners to restart state
real(kind=kind_real), pointer :: ud(:,:,:)
real(kind=kind_real), pointer :: vd(:,:,:)
real(kind=kind_real), pointer :: pe(:,:,:)
real(kind=kind_real), pointer :: pkz(:,:,:)
real(kind=kind_real), pointer :: pt(:,:,:)
real(kind=kind_real), pointer :: qils(:,:,:)
real(kind=kind_real), pointer :: qicn(:,:,:)
real(kind=kind_real), pointer :: qlls(:,:,:)
real(kind=kind_real), pointer :: qlcn(:,:,:)

! Pointers to background state
real(kind=kind_real), pointer :: ua(:,:,:)
real(kind=kind_real), pointer :: va(:,:,:)
real(kind=kind_real), pointer :: delp(:,:,:)
real(kind=kind_real), pointer :: ps(:,:,:)
real(kind=kind_real), pointer :: t(:,:,:)
real(kind=kind_real), pointer :: qi(:,:,:)
real(kind=kind_real), pointer :: ql(:,:,:)
real(kind=kind_real), pointer :: qilsf(:,:,:)
real(kind=kind_real), pointer :: qicnf(:,:,:)

real(kind=kind_real), allocatable :: pe_tmp(:,:,:)
real(kind=kind_real), target, allocatable :: pkz_tmp(:,:,:)

! Identity part of the change of variables
! ----------------------------------------
call copy_subset(xr%fields,xb%fields)


! D-Grid to A-Grid
! ----------------

if (self%do_wind) then

  call pointer_field_array(xr%fields, 'ud', ud)
  call pointer_field_array(xr%fields, 'vd', vd)

  call pointer_field_array(xb%fields, 'ua', ua)
  call pointer_field_array(xb%fields, 'va', va)

  call d2a(geom, ud, vd, ua, va)

endif

! Potential temperature to temperature
! ------------------------------------
if (self%do_temp) then

  call pointer_field_array(xr%fields, 'pt', pt)
  call pointer_field_array(xb%fields, 't' , t )

  if (.not. has_field(xr%fields,'pkz')) then
    if (has_field(xr%fields, 'delp')) then
      call pointer_field_array(xr%fields, 'delp' , delp )
      allocate(pe_tmp (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1))
      allocate(pkz_tmp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
      call delp_to_pe(geom, delp, pe_tmp)
      call pe_to_pk(geom, pe, pkz_tmp)
      pkz => pkz_tmp
    else
      call abor1_ftn("No way of getting pressures needed to convert temperature")
    endif
  else
    call pointer_field_array(xr%fields, 'pkz' , pkz )
  endif

  call pt_to_t(geom, pkz, pt, t)

endif

! Pressure to pressure thickness
! ------------------------------
if (self%do_pres) then

  call pointer_field_array(xr%fields, 'pe'  , pe)

  if (has_field(xb%fields, 'delp')) then
    call pointer_field_array(xb%fields, 'delp', delp)
    call pe_to_delp(geom,pe,delp)
  endif

  if (has_field(xb%fields, 'ps')) then
    call pointer_field_array(xb%fields, 'ps'  , ps)
    ps(:,:,1) = pe(:,:,geom%npz+1)
  endif

endif

! Four species of cloud to two species
! ------------------------------------

if (self%do_clds) then

  call pointer_field_array(xr%fields, 'qils', qils)
  call pointer_field_array(xr%fields, 'qicn', qicn)
  call pointer_field_array(xr%fields, 'qlls', qlls)
  call pointer_field_array(xr%fields, 'qlcn', qlcn)

  call pointer_field_array(xb%fields, 'ice_wat', qi)
  call pointer_field_array(xb%fields, 'liq_wat', ql)

  have_fractions = .true.
  if (.not.has_field(xb%fields, 'qilsf')) have_fractions = .false.
  if (.not.has_field(xb%fields, 'qicnf')) have_fractions = .false.

  if (have_fractions) then
    call pointer_field_array(xb%fields, 'qilsf', qilsf)
    call pointer_field_array(xb%fields, 'qicnf', qicnf)
  endif

  if (have_fractions) then
    call q4_to_q2(geom,qils,qicn,qlls,qlcn,qi,ql,qilsf,qicnf)
  else
    call q4_to_q2(geom,qils,qicn,qlls,qlcn,qi,ql)
  endif

endif

! Copy calendar infomation
! ------------------------
xb%calendar_type = xr%calendar_type
xb%date_init = xr%date_init

end subroutine changevar

! ------------------------------------------------------------------------------

subroutine changevarinverse(self,geom,xb,xr)

implicit none
type(fv3jedi_vc_geosrst2bkg), intent(inout) :: self
type(fv3jedi_geom),       intent(inout) :: geom
type(fv3jedi_state),      intent(in)    :: xb
type(fv3jedi_state),      intent(inout) :: xr

character(len=32) :: field_name

! Poitners to restart state
real(kind=kind_real), pointer :: ud(:,:,:)
real(kind=kind_real), pointer :: vd(:,:,:)
real(kind=kind_real), pointer :: pe(:,:,:)
real(kind=kind_real), pointer :: pkz(:,:,:)
real(kind=kind_real), pointer :: pt(:,:,:)
real(kind=kind_real), pointer :: qils(:,:,:)
real(kind=kind_real), pointer :: qicn(:,:,:)
real(kind=kind_real), pointer :: qlls(:,:,:)
real(kind=kind_real), pointer :: qlcn(:,:,:)

! Pointers to background state
real(kind=kind_real), pointer :: ua(:,:,:)
real(kind=kind_real), pointer :: va(:,:,:)
real(kind=kind_real), pointer :: delp(:,:,:)
real(kind=kind_real), pointer :: ps(:,:,:)
real(kind=kind_real), pointer :: t(:,:,:)
real(kind=kind_real), pointer :: qi(:,:,:)
real(kind=kind_real), pointer :: ql(:,:,:)
real(kind=kind_real), pointer :: qilsf(:,:,:)
real(kind=kind_real), pointer :: qicnf(:,:,:)

! Temporary arrays to hold pressures
real(kind=kind_real), allocatable :: pe_tmp(:,:,:)
real(kind=kind_real), allocatable :: pkz_tmp(:,:,:)


! Identity part of the change of variables
! ----------------------------------------
call copy_subset(xb%fields, xr%fields)

! A-Grid to D-Grid
! ----------------

if (self%do_wind) then

  call pointer_field_array(xr%fields, 'ud', ud)
  call pointer_field_array(xr%fields, 'vd', vd)

  call pointer_field_array(xb%fields, 'ua', ua)
  call pointer_field_array(xb%fields, 'va', va)

  call a2d(geom, ua, va, ud, vd)

endif

! Pressure to pressure thickness
! ------------------------------
if (self%do_pres) then

  call pointer_field_array(xr%fields, 'pe'  , pe )
  call pointer_field_array(xr%fields, 'pkz' , pkz )

  if (trim(self%pres_var) == 'delp') then

    call pointer_field_array(xb%fields, 'delp' , delp )
    call delp_to_pe(geom, delp, pe)

  elseif (trim(self%pres_var) == 'ps') then

    call pointer_field_array(xb%fields, 'ps' , ps )

    call ps_to_pe(geom, ps, pe)

  else

    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, must select a variable to set pe from. pres_var: ps or delp")

  endif

  ! Get p to the kappa
  call pe_to_pk(geom, pe, pkz)

endif

! Temperature to potential temperature
! ------------------------------------
if (self%do_temp) then

  call pointer_field_array(xr%fields, 'pt', pt)
  call pointer_field_array(xb%fields, 't' , t )

  if (.not. self%do_pres) then

    allocate(pe_tmp (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1))
    allocate(pkz_tmp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))

    call pointer_field_array(xb%fields, 'delp' , delp )

    call delp_to_pe(geom, delp, pe_tmp)
    call pe_to_pk(geom, pe_tmp, pkz_tmp)

    call t_to_pt(geom, pkz_tmp, t, pt)

    deallocate(pe_tmp,pkz_tmp)

  else

    call t_to_pt(geom, pkz, t, pt)

  endif

endif

! Four species of cloud to two species
! ------------------------------------

if (self%do_clds) then

  call pointer_field_array(xr%fields, 'qils', qils)
  call pointer_field_array(xr%fields, 'qicn', qicn)
  call pointer_field_array(xr%fields, 'qlls', qlls)
  call pointer_field_array(xr%fields, 'qlcn', qlcn)

  call pointer_field_array(xb%fields, 'ice_wat', qi)
  call pointer_field_array(xb%fields, 'liq_wat', ql)
  call pointer_field_array(xb%fields, 'qilsf', qilsf)
  call pointer_field_array(xb%fields, 'qicnf', qicnf)

  call q2_to_q4(geom, qi, ql, qilsf, qicnf, qils, qicn, qlls, qlcn)

endif

! Copy calendar infomation
! ------------------------
xr%calendar_type = xb%calendar_type
xr%date_init = xb%date_init

end subroutine changevarinverse

! ------------------------------------------------------------------------------

end module fv3jedi_vc_geosrst2bkg_mod
