! (C) Copyright 2018-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_varcha_a2m_mod

! iso
use iso_c_binding

! fckit
use fckit_configuration_module, only: fckit_configuration

! oops
use datetime_mod

! fv3-jedi
use fv3jedi_kinds_mod,                  only: kind_real
use fv3jedi_geom_mod,                   only: fv3jedi_geom
use fv3jedi_state_mod,                  only: fv3jedi_state
use fv3jedi_field_mod,                  only: copy_subset
use wind_vt_mod,                        only: a_to_d, d_to_a

implicit none
private
public :: fv3jedi_varcha_a2m
type :: fv3jedi_varcha_a2m
  contains
    procedure :: create
    procedure :: delete
    procedure :: changevar
    procedure :: changevarinverse
end type fv3jedi_varcha_a2m

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, conf)

class(fv3jedi_varcha_a2m), intent(inout) :: self
type(fv3jedi_geom),        intent(inout) :: geom
type(fckit_configuration), intent(in)    :: conf

! Nothing to do

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

class(fv3jedi_varcha_a2m), intent(inout) :: self

! Nothing to do

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine changevar(self,geom,xana,xmod)

class(fv3jedi_varcha_a2m), intent(inout) :: self
type(fv3jedi_geom),        intent(inout) :: geom
type(fv3jedi_state),       intent(in)    :: xana
type(fv3jedi_state),       intent(inout) :: xmod

integer :: index_ana, index_mod, index_ana_found
integer :: k
logical :: failed

real(kind=kind_real), pointer, dimension(:,:,:) :: xana_ua
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_va

real(kind=kind_real), pointer, dimension(:,:,:) :: xmod_ud
real(kind=kind_real), pointer, dimension(:,:,:) :: xmod_vd

do index_mod = 1, xmod%nf

  index_ana_found = -1
  failed = .true.

  !Check analysis for presence of field
  do index_ana = 1, xana%nf
    if (xmod%fields(index_mod)%short_name == xana%fields(index_ana)%short_name) then
      index_ana_found = index_ana
      exit
    endif
  enddo

  if (index_ana_found >= 0) then

    !OK, direct copy
    xmod%fields(index_mod)%array = xana%fields(index_ana_found)%array
    failed = .false.

  elseif (xmod%fields(index_mod)%short_name == 'ud') then

    !Special case: A-grid analysis, D-Grid model
    if (xana%has_field('ua')) then
      call xana%get_field('ua', xana_ua)
      call xana%get_field('va', xana_va)
      call xmod%get_field('ud', xmod_ud)
      call xmod%get_field('vd', xmod_vd)
      call a_to_d(geom, xana_ua, xana_va, xmod_ud, xmod_vd)
      xmod_ud(:,geom%jec+1,:) = 0.0_kind_real
      xmod_vd(geom%iec+1,:,:) = 0.0_kind_real
      failed = .false.
    endif

  elseif (xmod%fields(index_mod)%short_name == 'vd') then

    !Special case: A-grid analysis, D-Grid model
    if (xana%has_field('ua')) then
      !Already done above
      failed = .false.
    endif

  endif

  if (failed) &
    call abor1_ftn("fv3jedi_varcha_a2m_mod.changevar: Found no way of getting " &
                   // trim(xmod%fields(index_mod)%short_name) // " from the analysis state.")

enddo

end subroutine changevar

! --------------------------------------------------------------------------------------------------

subroutine changevarinverse(self,geom,xmod,xana)

class(fv3jedi_varcha_a2m), intent(inout) :: self
type(fv3jedi_geom),        intent(inout) :: geom
type(fv3jedi_state),       intent(in)    :: xmod
type(fv3jedi_state),       intent(inout) :: xana

integer :: index_ana, index_mod, index_mod_found
logical :: failed

real(kind=kind_real), pointer, dimension(:,:,:) :: xana_ua
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_va

real(kind=kind_real), pointer, dimension(:,:,:) :: xmod_ud
real(kind=kind_real), pointer, dimension(:,:,:) :: xmod_vd

do index_ana = 1, xana%nf

  index_mod_found = -1
  failed = .true.

  !Check model for presence of field
  do index_mod = 1, xmod%nf
    if (xana%fields(index_ana)%short_name == xmod%fields(index_mod)%short_name) then
      index_mod_found = index_mod
      exit
    endif
  enddo

  if (index_mod_found >= 0) then

    !OK, direct copy
    failed = .false.
    xana%fields(index_ana)%array = xmod%fields(index_mod_found)%array

  elseif (xana%fields(index_ana)%short_name == 'ua') then

    !Special case: A-grid analysis, D-Grid model
    if (xmod%has_field('ud')) then
      call xana%get_field('ua', xana_ua)
      call xana%get_field('va', xana_va)
      call xmod%get_field('ud', xmod_ud)
      call xmod%get_field('vd', xmod_vd)
      xmod_ud(:,geom%jec+1,:) = 0.0_kind_real
      xmod_vd(geom%iec+1,:,:) = 0.0_kind_real
      call d_to_a(geom, xmod_ud, xmod_vd, xana_ua, xana_va)
      failed = .false.
    endif

  elseif (xana%fields(index_ana)%short_name == 'va') then

    !Special case: A-grid analysis, D-Grid model
    if (xmod%has_field('ud')) then
      !Already done above
      failed = .false.
    endif

  endif

  if (failed) &
    call abor1_ftn("fv3jedi_varcha_a2m_mod.changevarinverse: Found no way of getting " &
                   // trim(xmod%fields(index_mod)%short_name) // " from the model state.")

enddo

end subroutine changevarinverse

! --------------------------------------------------------------------------------------------------

end module fv3jedi_varcha_a2m_mod
