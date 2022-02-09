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
use fv3jedi_io_fms_mod,                 only: fv3jedi_io_fms
use fv3jedi_io_cube_sphere_history_mod, only: fv3jedi_io_cube_sphere_history
use fv3jedi_field_mod,                  only: copy_subset
use wind_vt_mod,                        only: a_to_d, d_to_a

implicit none
private
public :: fv3jedi_varcha_a2m
type :: fv3jedi_varcha_a2m
  character(len=96)     :: filetype !IO type
  type(fv3jedi_io_fms)  :: fms
  type(fv3jedi_io_cube_sphere_history) :: csh
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

implicit none
class(fv3jedi_varcha_a2m), intent(inout) :: self
type(fv3jedi_geom),        intent(inout) :: geom
type(fckit_configuration), intent(in)    :: conf

character(len=:), allocatable :: str

! Get IO type to use
call conf%get_or_die("filetype", str)
self%filetype = str

! Setup IO from config
if (trim(self%filetype) == 'fms restart') then
  call self%fms%create(conf, geom%domain, geom%npz)
elseif (trim(self%filetype) == 'cube sphere history') then
  call self%csh%create(geom, conf)
else
  call abor1_ftn("fv3jedi_varcha_a2m_mod: filetype must be cube sphere history or fms restart, is "// &
                 trim(self%filetype)//".")
endif

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)

implicit none
class(fv3jedi_varcha_a2m), intent(inout) :: self

if (trim(self%filetype) == 'fms restart') call self%fms%delete()
if (trim(self%filetype) == 'cube sphere history') call self%csh%delete()

end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine changevar(self,geom,xana,xmod)

implicit none
class(fv3jedi_varcha_a2m), intent(inout) :: self
type(fv3jedi_geom),        intent(inout) :: geom
type(fv3jedi_state),       intent(in)    :: xana
type(fv3jedi_state),       intent(inout) :: xmod

integer :: index_ana, index_mod, index_ana_found
integer :: k
logical :: failed

real(kind=kind_real), pointer, dimension(:,:,:) :: xana_ua
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_va
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_ps

real(kind=kind_real), pointer, dimension(:,:,:) :: xmod_ud
real(kind=kind_real), pointer, dimension(:,:,:) :: xmod_vd
real(kind=kind_real), pointer, dimension(:,:,:) :: xmod_delp

do index_mod = 1, xmod%nf

  index_ana_found = -1
  failed = .true.

  !Check analysis for presence of field
  do index_ana = 1, xana%nf
    if (xmod%fields(index_mod)%fv3jedi_name == xana%fields(index_ana)%fv3jedi_name) then
      index_ana_found = index_ana
      exit
    endif
  enddo

  if (index_ana_found >= 0) then

    !OK, direct copy
    xmod%fields(index_mod)%array = xana%fields(index_ana_found)%array
    failed = .false.
    if (xmod%f_comm%rank() == 0) write(*,"(A, A10, A, A10)") &
        "A2M ChangeVar: analysis state "//xana%fields(index_ana_found)%fv3jedi_name(1:10)&
        //" => model state "//xmod%fields(index_mod)%fv3jedi_name(1:10)

  elseif (xmod%fields(index_mod)%fv3jedi_name == 'ud') then

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
      if (xmod%f_comm%rank() == 0) write(*,"(A)") &
          "A2M ChangeVar: analysis state ua         => model state ud"
    endif

  elseif (xmod%fields(index_mod)%fv3jedi_name == 'vd') then

    !Special case: A-grid analysis, D-Grid model
    if (xana%has_field('ua')) then
      !Already done above
      failed = .false.
      if (xmod%f_comm%rank() == 0) write(*,"(A)") &
          "A2M ChangeVar: analysis state va         => model state vd"
    endif

  elseif (xmod%fields(index_mod)%fv3jedi_name == 'delp') then

    !Special case: ps in analysis, delp in model
    if (xana%has_field('ps')) then
      call xana%get_field('ps',   xana_ps)
      call xmod%get_field('delp', xmod_delp)
      do k = 1,geom%npz
        xmod_delp(:,:,k) = (geom%ak(k+1)-geom%ak(k)) + (geom%bk(k+1)-geom%bk(k))*xana_ps(:,:,1)
      enddo
      failed = .false.
      if (xmod%f_comm%rank() == 0) write(*,"(A)") &
          "A2M ChangeVar: analysis state ps         => model state delp"
    endif

  endif

  if (failed) then

    if (xmod%f_comm%rank() == 0) write(*,"(A)") &
        "Found no way of getting "//trim(xmod%fields(index_mod)%fv3jedi_name)//" from the analysis state."//&
        "Attempting to read from file"

    if (trim(self%filetype) == 'fms restart') then
      call self%fms%read(xmod%time, xmod%fields(index_mod:index_mod))
    elseif (trim(self%filetype) == 'cube sphere history') then
      call self%csh%read(xmod%time, xmod%fields(index_mod:index_mod))
    endif

  endif

enddo

end subroutine changevar

! --------------------------------------------------------------------------------------------------

subroutine changevarinverse(self,geom,xmod,xana)

implicit none
class(fv3jedi_varcha_a2m), intent(inout) :: self
type(fv3jedi_geom),        intent(inout) :: geom
type(fv3jedi_state),       intent(in)    :: xmod
type(fv3jedi_state),       intent(inout) :: xana

integer :: index_ana, index_mod, index_mod_found
logical :: failed

real(kind=kind_real), pointer, dimension(:,:,:) :: xana_ua
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_va
real(kind=kind_real), pointer, dimension(:,:,:) :: xana_ps

real(kind=kind_real), pointer, dimension(:,:,:) :: xmod_ud
real(kind=kind_real), pointer, dimension(:,:,:) :: xmod_vd
real(kind=kind_real), pointer, dimension(:,:,:) :: xmod_delp

do index_ana = 1, xana%nf

  index_mod_found = -1
  failed = .true.

  !Check analysis for presence of field
  do index_mod = 1, xmod%nf
    if (xana%fields(index_ana)%fv3jedi_name == xmod%fields(index_mod)%fv3jedi_name) then
      index_mod_found = index_mod
      exit
    endif
  enddo

  if (index_mod_found >= 0) then

    !OK, direct copy
    failed = .false.
    xana%fields(index_ana)%array = xmod%fields(index_mod_found)%array
    if (xana%f_comm%rank() == 0) write(*,"(A, A10, A, A10)") &
        "A2M ChangeVarInverse: model state "//xmod%fields(index_mod_found)%fv3jedi_name(1:10)&
        //" => analysis state "//xana%fields(index_ana)%fv3jedi_name(1:10)

  elseif (xana%fields(index_ana)%fv3jedi_name == 'ua') then

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
      if (xana%f_comm%rank() == 0) write(*,"(A)") &
          "A2M ChangeVarInverse: model state ud         => analysis state ua"
    endif

  elseif (xana%fields(index_ana)%fv3jedi_name == 'va') then

    !Special case: A-grid analysis, D-Grid model
    if (xmod%has_field('ud')) then
      !Already done above
      failed = .false.
      if (xana%f_comm%rank() == 0) write(*,"(A)") &
          "A2M ChangeVarInverse: model state vd         => analysis state va"
    endif

  elseif (xana%fields(index_ana)%fv3jedi_name == 'ps') then

    !Special case: ps in analysis, delp in model
    if (xmod%has_field('delp')) then
      call xana%get_field('ps',   xana_ps)
      call xmod%get_field('delp', xmod_delp)
      xana_ps(:,:,1) = sum(xmod_delp,3)
      failed = .false.
      if (xana%f_comm%rank() == 0) write(*,"(A)") &
          "A2M ChangeVarInverse: model state delp       => analysis state ps"
    endif

  endif

  if (failed) call abor1_ftn("fv3jedi_linvarcha_a2m_mod.changevarinverse: found no way of getting "//&
                             trim(xana%fields(index_ana)%fv3jedi_name)//" from the model state" )

enddo

end subroutine changevarinverse

! --------------------------------------------------------------------------------------------------

end module fv3jedi_varcha_a2m_mod
