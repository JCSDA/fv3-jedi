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
use fv3jedi_io_gfs_mod,  only: fv3jedi_io_gfs
use fv3jedi_io_geos_mod, only: fv3jedi_io_geos

use fv3jedi_field_mod, only: copy_subset, get_field_array

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

character(len=32) :: field_name
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

logical :: have_pkz
real(kind=kind_real), allocatable :: pe_tmp(:,:,:)
real(kind=kind_real), allocatable :: pkz_tmp(:,:,:)

! Identity part of the change of variables
! ----------------------------------------
call copy_subset(xr%fields,xb%fields)


! D-Grid to A-Grid
! ----------------

if (self%do_wind) then

  field_name = 'ud'
  if (.not. get_field_array(xr%fields, trim(field_name), ud )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevar, no "//trim(field_name)//" in inputs")
  field_name = 'vd'
  if (.not. get_field_array(xr%fields, trim(field_name), vd )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevar, no "//trim(field_name)//" in inputs")

  field_name = 'ua'
  if (.not. get_field_array(xb%fields, trim(field_name), ua )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevar, no "//trim(field_name)//" in outputs")
  field_name = 'va'
  if (.not. get_field_array(xb%fields, trim(field_name), va )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevar, no "//trim(field_name)//" in outputs")

  call d2a(geom, ud, vd, ua, va)

endif

! Potential temperature to temperature
! ------------------------------------
if (self%do_temp) then

  field_name = 'pt'
  if (.not. get_field_array(xr%fields, trim(field_name), pt )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevar, no "//trim(field_name)//" in inputs")

  field_name = 't'
  if (.not. get_field_array(xb%fields, trim(field_name), t )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevar, no "//trim(field_name)//" in outputs")

  have_pkz = get_field_array(xr%fields, 'pkz', pkz )

  if (.not. have_pkz) then
    if (get_field_array(xr%fields, 'delp', delp )) then
      allocate(pe_tmp (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1))
      allocate(pkz_tmp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
      call delp_to_pe(geom, delp, pe_tmp)
      call pe_to_pk(geom, pe, pkz_tmp)
      pkz = pkz_tmp
    else
      call abor1_ftn("No way of getting pressures needed to convert temperature")
    endif
  endif

  call pt_to_t(geom,pkz,pt,t)

endif

! Pressure to pressure thickness
! ------------------------------
if (self%do_pres) then

  field_name = 'pe'
  if (.not. get_field_array(xr%fields, trim(field_name), pe )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevar, no "//trim(field_name)//" in inputs")

  field_name = 'delp'
  if (get_field_array(xb%fields, trim(field_name), delp )) then
    call pe_to_delp(geom,pe,delp)
  endif

  field_name = 'ps'
  if (get_field_array(xb%fields, trim(field_name), ps )) then
    ps(:,:,1) = pe(:,:,geom%npz+1)
  endif

endif

! Four species of cloud to two species
! ------------------------------------

if (self%do_clds) then

  field_name = 'qils'
  if (.not. get_field_array(xr%fields, trim(field_name), qils )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevar, no "//trim(field_name)//" in inputs")
  field_name = 'qicn'
  if (.not. get_field_array(xr%fields, trim(field_name), qicn )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevar, no "//trim(field_name)//" in inputs")
  field_name = 'qlls'
  if (.not. get_field_array(xr%fields, trim(field_name), qlls )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevar, no "//trim(field_name)//" in inputs")
  field_name = 'qlcn'
  if (.not. get_field_array(xr%fields, trim(field_name), qlcn )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevar, no "//trim(field_name)//" in inputs")

  field_name = 'qi'
  if (.not. get_field_array(xb%fields, trim(field_name), qi )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevar, no "//trim(field_name)//" in outputs")
  field_name = 'ql'
  if (.not. get_field_array(xb%fields, trim(field_name), ql )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevar, no "//trim(field_name)//" in outputs")

  have_fractions = .true.
  field_name = 'qilsf'
  if (.not. get_field_array(xb%fields, trim(field_name), qilsf )) have_fractions = .false.
  field_name = 'qicnf'
  if (.not. get_field_array(xb%fields, trim(field_name), qicnf )) have_fractions = .false.

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
call copy_subset(xb%fields,xr%fields)

! A-Grid to D-Grid
! ----------------

if (self%do_wind) then

  field_name = 'ud'
  if (.not. get_field_array(xr%fields, trim(field_name), ud )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in outputs")
  field_name = 'vd'
  if (.not. get_field_array(xr%fields, trim(field_name), vd )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in outputs")

  field_name = 'ua'
  if (.not. get_field_array(xb%fields, trim(field_name), ua )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in inputs")
  field_name = 'va'
  if (.not. get_field_array(xb%fields, trim(field_name), va )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in inputs")

  call a2d(geom, ua, va, ud, vd)

endif

! Pressure to pressure thickness
! ------------------------------
if (self%do_pres) then

  field_name = 'pe'
  if (.not. get_field_array(xr%fields, trim(field_name), pe )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in outputs")
  field_name = 'pkz'
  if (.not. get_field_array(xr%fields, trim(field_name), pkz )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in outputs")

  if (trim(self%pres_var) == 'delp') then

    field_name = 'delp'
    if (.not. get_field_array(xb%fields, trim(field_name), delp )) &
      call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in inputs")

    call delp_to_pe(geom, delp, pe)

  elseif (trim(self%pres_var) == 'ps') then

    field_name = 'ps'
    if (.not. get_field_array(xb%fields, trim(field_name), ps )) &
      call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in inputs")

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

  field_name = 'pt'
  if (.not. get_field_array(xr%fields, trim(field_name), pt )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in outputs")

  field_name = 't'
  if (.not. get_field_array(xb%fields, trim(field_name), t )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in inputs")

  if (.not. self%do_pres) then

    allocate(pe_tmp (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1))
    allocate(pkz_tmp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))

    field_name = 'delp'
    if (.not. get_field_array(xb%fields, trim(field_name), delp )) &
      call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in inputs")

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

  field_name = 'qils'
  if (.not. get_field_array(xr%fields, trim(field_name), qils )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in outputs")
  field_name = 'qicn'
  if (.not. get_field_array(xr%fields, trim(field_name), qicn )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in outputs")
  field_name = 'qlls'
  if (.not. get_field_array(xr%fields, trim(field_name), qlls )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in outputs")
  field_name = 'qlcn'
  if (.not. get_field_array(xr%fields, trim(field_name), qlcn )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in outputs")

  field_name = 'qi'
  if (.not. get_field_array(xb%fields, trim(field_name), qi )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in inputs")
  field_name = 'ql'
  if (.not. get_field_array(xb%fields, trim(field_name), ql )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in inputs")
  field_name = 'qilsf'
  if (.not. get_field_array(xb%fields, trim(field_name), qilsf )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in inputs")
  field_name = 'qicnf'
  if (.not. get_field_array(xb%fields, trim(field_name), qicnf )) &
    call abor1_ftn("fv3jedi_vc_geosrst2bkg_mod.changevarinverse, no "//trim(field_name)//" in inputs")

  call q2_to_q4(geom,qi,ql,qilsf,qicnf,qils,qicn,qlls,qlcn)

endif


! Copy calendar infomation
! ------------------------
xr%calendar_type = xb%calendar_type
xr%date_init = xb%date_init

end subroutine changevarinverse

! ------------------------------------------------------------------------------

end module fv3jedi_vc_geosrst2bkg_mod
