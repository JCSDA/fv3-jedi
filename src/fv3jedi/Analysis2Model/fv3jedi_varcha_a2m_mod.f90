! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_varcha_a2m_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use datetime_mod

use fckit_log_module, only : fckit_log

use fv3jedi_kinds_mod,   only: kind_real
use fv3jedi_geom_mod,    only: fv3jedi_geom
use fv3jedi_state_mod,   only: fv3jedi_state
use fv3jedi_io_gfs_mod,  only: fv3jedi_io_gfs
use fv3jedi_io_geos_mod, only: fv3jedi_io_geos

use wind_vt_mod, only: a2d, d2a

implicit none
private

public :: fv3jedi_varcha_a2m
public :: create
public :: delete
public :: changevar
public :: changevarinverse

type :: fv3jedi_varcha_a2m
  character(len=10)  :: filetype !Model type
  type(fv3jedi_io_gfs)  :: gfs
  type(fckit_configuration) :: f_conf
end type fv3jedi_varcha_a2m

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, c_conf)

implicit none
type(fv3jedi_varcha_a2m), intent(inout) :: self
type(fv3jedi_geom),       intent(inout) :: geom
type(c_ptr),              intent(in)    :: c_conf

character(len=:), allocatable :: str

! Fortran configuration
self%f_conf = fckit_configuration(c_conf)

! Model type
call self%f_conf%get_or_die("filetype",str)
self%filetype = str
deallocate(str)

if (trim(self%filetype) == 'gfs') then
  call self%gfs%setup(self%f_conf)
endif

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_varcha_a2m), intent(inout) :: self

end subroutine delete

! ------------------------------------------------------------------------------

subroutine changevar(self,geom,xana,xmod,vdt)

implicit none
type(fv3jedi_varcha_a2m), intent(inout) :: self
type(fv3jedi_geom),       intent(inout) :: geom
type(fv3jedi_state),      intent(in)    :: xana
type(fv3jedi_state),      intent(inout) :: xmod
type(datetime),           intent(inout) :: vdt

integer :: index_ana, index_mod, index_ana_found
integer :: k
logical :: failed

type(fv3jedi_io_geos) :: geos

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
        "A2M ChangeVar: analysis state "//xana%fields(index_ana_found)%fv3jedi_name&
        //" => model state "//xmod%fields(index_mod)%fv3jedi_name

  elseif (xmod%fields(index_mod)%fv3jedi_name == 'ud') then

    !Special case: A-grid analysis, D-Grid model
    if (associated(xana%ua) .and. associated(xana%va)) then
      call a2d(geom, xana%ua, xana%va, xmod%ud, xmod%vd)
      xmod%ud(:,geom%jec+1,:) = 0.0_kind_real
      xmod%vd(geom%iec+1,:,:) = 0.0_kind_real
      failed = .false.
      if (xmod%f_comm%rank() == 0) write(*,"(A)") &
          "A2M ChangeVar: analysis state ua         => model state ud"
    endif

  elseif (xmod%fields(index_mod)%fv3jedi_name == 'vd') then

    !Special case: A-grid analysis, D-Grid model
    if (associated(xana%ua) .and. associated(xana%va)) then
      !Already done above
      failed = .false.
      if (xmod%f_comm%rank() == 0) write(*,"(A)") &
          "A2M ChangeVar: analysis state va         => model state vd"
    endif

  elseif (xmod%fields(index_mod)%fv3jedi_name == 'delp') then

    !Special case: ps in analysis, delp in model
    if (associated(xana%ps)) then
      do k = 1,geom%npz
        xmod%delp(:,:,k) = (geom%ak(k+1)-geom%ak(k)) + (geom%bk(k+1)-geom%bk(k))*xana%ps(:,:,1)
      enddo
      failed = .false.
      if (xmod%f_comm%rank() == 0) write(*,"(A)") &
          "A2M ChangeVar: analysis state ps         => model state delp"
    endif

  endif

  if (failed) then

    if (xmod%f_comm%rank() == 0) write(*,"(A)") &
        "Found no way of getting "//xmod%fields(index_mod)%fv3jedi_name//" from the analysis state."//&
        "Attempting to read from file"

    if (trim(self%filetype) == 'gfs') then
      call self%gfs%read_fields( geom, xmod%fields(index_mod:index_mod) )
    elseif (trim(self%filetype) == 'geos') then
      call geos%setup(geom, xmod%fields(index_mod:index_mod), vdt, 'read', self%f_conf)
      call geos%read_fields(geom, xmod%fields(index_mod:index_mod))
      call geos%delete()
    endif

  endif

enddo

end subroutine changevar

! ------------------------------------------------------------------------------

subroutine changevarinverse(self,geom,xmod,xana,vdt)

implicit none
type(fv3jedi_varcha_a2m), intent(inout) :: self
type(fv3jedi_geom),       intent(inout) :: geom
type(fv3jedi_state),      intent(in)    :: xmod
type(fv3jedi_state),      intent(inout) :: xana
type(datetime),           intent(inout) :: vdt

integer :: index_ana, index_mod, index_mod_found
logical :: failed

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
        "A2M ChangeVarInverse: model state "//xmod%fields(index_mod_found)%fv3jedi_name&
        //" => analysis state "//xana%fields(index_ana)%fv3jedi_name

  elseif (xana%fields(index_ana)%fv3jedi_name == 'ua') then

    !Special case: A-grid analysis, D-Grid model
    if (associated(xmod%ud) .and. associated(xmod%vd)) then
      xmod%ud(:,geom%jec+1,:) = 0.0_kind_real
      xmod%vd(geom%iec+1,:,:) = 0.0_kind_real
      call d2a(geom, xmod%ud, xmod%vd, xana%ua, xana%va)
      failed = .false.
      if (xana%f_comm%rank() == 0) write(*,"(A)") &
          "A2M ChangeVarInverse: model state ud         => analysis state ua"
    endif

  elseif (xana%fields(index_ana)%fv3jedi_name == 'va') then

    !Special case: A-grid analysis, D-Grid model
    if (associated(xmod%ud) .and. associated(xmod%vd)) then
      !Already done above
      failed = .false.
      if (xana%f_comm%rank() == 0) write(*,"(A)") &
          "A2M ChangeVarInverse: model state vd         => analysis state va"
    endif

  elseif (xana%fields(index_ana)%fv3jedi_name == 'ps') then

    !Special case: ps in analysis, delp in model
    if (associated(xmod%delp)) then
      xana%ps(:,:,1) = sum(xmod%delp,3)
      failed = .false.
      if (xana%f_comm%rank() == 0) write(*,"(A)") &
          "A2M ChangeVarInverse: model state delp       => analysis state ps"
    endif

  endif

  if (failed) call abor1_ftn("fv3jedi_linvarcha_a2m_mod.changevarinverse: found no way of getting "//&
                             xana%fields(index_ana)%fv3jedi_name//" from the model state" )

enddo


end subroutine changevarinverse

! ------------------------------------------------------------------------------

end module fv3jedi_varcha_a2m_mod
