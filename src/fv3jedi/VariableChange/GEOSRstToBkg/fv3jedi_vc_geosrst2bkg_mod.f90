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

use fv3jedi_field_mod, only: copy_subset

use wind_vt_mod, only: a_to_d, d_to_a
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

if( .not. conf%get('do_wind',        self%do_wind) ) self%do_wind = .true.

if( .not. conf%get('do_temperature', self%do_temp) ) self%do_temp = .true.

if( .not. conf%get('do_pressure',    self%do_pres) ) self%do_pres = .true.

if( .not. conf%get('do_clouds',      self%do_clds) ) self%do_clds = .true.

self%pres_var = 'delp'
if( conf%has('pres_var') ) then
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

  call xr%get_field('ud', ud)
  call xr%get_field('vd', vd)

  call xb%get_field('ua', ua)
  call xb%get_field('va', va)

  call d_to_a(geom, ud, vd, ua, va)

endif

! Potential temperature to temperature
! ------------------------------------
if (self%do_temp) then

  call xr%get_field('pt', pt)
  call xb%get_field('t' , t )

  if (.not. xr%has_field('pkz')) then
    if (xr%has_field( 'delp')) then
      call xr%get_field('delp' , delp )
      allocate(pe_tmp (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1))
      allocate(pkz_tmp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))
      call delp_to_pe(geom, delp, pe_tmp)
      call pe_to_pk(geom, pe, pkz_tmp)
      pkz => pkz_tmp
    else
      call abor1_ftn("No way of getting pressures needed to convert temperature")
    endif
  else
    call xr%get_field('pkz' , pkz )
  endif

  call pt_to_t(geom, pkz, pt, t)

endif

! Pressure to pressure thickness
! ------------------------------
if (self%do_pres) then

  call xr%get_field('pe'  , pe)

  if (xb%has_field( 'delp')) then
    call xb%get_field('delp', delp)
    call pe_to_delp(geom,pe,delp)
  endif

  if (xb%has_field( 'ps')) then
    call xb%get_field('ps'  , ps)
    ps(:,:,1) = pe(:,:,geom%npz+1)
  endif

endif

! Four species of cloud to two species
! ------------------------------------

if (self%do_clds) then

  call xr%get_field('qils', qils)
  call xr%get_field('qicn', qicn)
  call xr%get_field('qlls', qlls)
  call xr%get_field('qlcn', qlcn)

  call xb%get_field('ice_wat', qi)
  call xb%get_field('liq_wat', ql)

  have_fractions = .true.
  if (.not.xb%has_field( 'qilsf')) have_fractions = .false.
  if (.not.xb%has_field( 'qicnf')) have_fractions = .false.

  if (have_fractions) then
    call xb%get_field('qilsf', qilsf)
    call xb%get_field('qicnf', qicnf)
  endif

  if (have_fractions) then
    call q4_to_q2(geom,qils,qicn,qlls,qlcn,qi,ql,qilsf,qicnf)
  else
    call q4_to_q2(geom,qils,qicn,qlls,qlcn,qi,ql)
  endif

endif

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

  call xr%get_field('ud', ud)
  call xr%get_field('vd', vd)

  call xb%get_field('ua', ua)
  call xb%get_field('va', va)

  call a_to_d(geom, ua, va, ud, vd)

endif

! Pressure to pressure thickness
! ------------------------------
if (self%do_pres) then

  call xr%get_field('pe'  , pe )
  call xr%get_field('pkz' , pkz )

  if (trim(self%pres_var) == 'delp') then

    call xb%get_field('delp' , delp )
    call delp_to_pe(geom, delp, pe)

  elseif (trim(self%pres_var) == 'ps') then

    call xb%get_field('ps' , ps )

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

  call xr%get_field('pt', pt)
  call xb%get_field('t' , t )

  if (.not. self%do_pres) then

    allocate(pe_tmp (geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz+1))
    allocate(pkz_tmp(geom%isc:geom%iec,geom%jsc:geom%jec,1:geom%npz))

    call xb%get_field('delp' , delp )

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

  call xr%get_field('qils', qils)
  call xr%get_field('qicn', qicn)
  call xr%get_field('qlls', qlls)
  call xr%get_field('qlcn', qlcn)

  call xb%get_field('ice_wat', qi)
  call xb%get_field('liq_wat', ql)
  call xb%get_field('qilsf', qilsf)
  call xb%get_field('qicnf', qicnf)

  call q2_to_q4(geom, qi, ql, qilsf, qicnf, qils, qicn, qlls, qlcn)

endif

end subroutine changevarinverse

! ------------------------------------------------------------------------------

end module fv3jedi_vc_geosrst2bkg_mod
