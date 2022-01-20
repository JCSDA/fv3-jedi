! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_linvarcha_a2m_mod

use iso_c_binding
use fckit_configuration_module, only: fckit_configuration
use fv3jedi_kinds_mod

use fv3jedi_geom_mod,      only: fv3jedi_geom
use fv3jedi_state_mod,     only: fv3jedi_state
use fv3jedi_increment_mod, only: fv3jedi_increment

use fv3jedi_field_mod, only: copy_subset

use wind_vt_mod, only: a_to_d, d_to_a, d_to_a_ad, a_to_d_ad

implicit none
private

public :: fv3jedi_linvarcha_a2m
public :: create
public :: delete
public :: multiply
public :: multiplyadjoint
public :: multiplyinverse
public :: multiplyinverseadjoint

type :: fv3jedi_linvarcha_a2m
 integer :: dummy
end type fv3jedi_linvarcha_a2m

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine create(self, geom, bg, fg, c_conf)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_state),         intent(in)    :: bg
type(fv3jedi_state),         intent(in)    :: fg
type(c_ptr),                 intent(in)    :: c_conf

type(fckit_configuration) :: f_conf

f_conf = fckit_configuration(c_conf)

!Nothing yet as transforms are linear

end subroutine create

! ------------------------------------------------------------------------------

subroutine delete(self)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self

end subroutine delete

! ------------------------------------------------------------------------------

subroutine multiply(self,geom,xana,xmod)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: xana
type(fv3jedi_increment),     intent(inout) :: xmod

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
        "A2M Multiply: analysis increment "//xana%fields(index_ana_found)%fv3jedi_name(1:10)&
        //" => linearized model "//xmod%fields(index_mod)%fv3jedi_name(1:10)

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
          "A2M Multiply: analysis increment ua         => linearized model ud"
    endif

  elseif (xmod%fields(index_mod)%fv3jedi_name == 'vd') then

    !Special case: A-grid analysis, D-Grid model
    if (xana%has_field('ua')) then
      !Already done above
      failed = .false.
      if (xmod%f_comm%rank() == 0) write(*,"(A)") &
          "A2M Multiply: analysis increment va         => linearized model ud"
    endif

  elseif (xmod%fields(index_mod)%fv3jedi_name == 'delp') then

    !Special case: ps in analysis, delp in model
    if (xana%has_field('ps')) then
      call xana%get_field('ps',   xana_ps)
      call xmod%get_field('delp', xmod_delp)
      do k = 1,geom%npz
        xmod_delp(:,:,k) = (geom%bk(k+1)-geom%bk(k))*xana_ps(:,:,1)
      enddo
      failed = .false.
      if (xmod%f_comm%rank() == 0) write(*,"(A)") &
          "A2M Multiply: analysis increment ps         => linearized model ud"
    endif

  endif

  if (failed) call abor1_ftn("fv3jedi_linvarcha_a2m_mod.multiply: found no way of getting "//&
                             trim(xmod%fields(index_mod)%fv3jedi_name)//" from the analysis increment" )

enddo

end subroutine multiply

! ------------------------------------------------------------------------------

subroutine multiplyadjoint(self,geom,xmod,xana)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: xmod
type(fv3jedi_increment),     intent(inout) :: xana

integer :: index_ana, index_mod, index_mod_found
integer :: k
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
    xana%fields(index_ana)%array = xmod%fields(index_mod_found)%array
    failed = .false.
    if (xana%f_comm%rank() == 0) write(*,"(A, A10, A, A10)") &
        "A2M MultiplyAdjoint: linearized model "//xmod%fields(index_mod_found)%fv3jedi_name(1:10)&
        //" => analysis increment "//xana%fields(index_ana)%fv3jedi_name(1:10)

  elseif (xana%fields(index_ana)%fv3jedi_name == 'ua') then

    !Special case: A-grid analysis, D-Grid model
    if (xmod%has_field('ud')) then
      call xana%get_field('ua', xana_ua)
      call xana%get_field('va', xana_va)
      call xmod%get_field('ud', xmod_ud)
      call xmod%get_field('vd', xmod_vd)
      xmod_ud(:,geom%jec+1,:) = 0.0_kind_real
      xmod_vd(geom%iec+1,:,:) = 0.0_kind_real
      call a_to_d_ad(geom, xana_ua, xana_va, xmod_ud, xmod_vd)
      failed = .false.
      if (xana%f_comm%rank() == 0) write(*,"(A)") &
          "A2M MultiplyAdjoint: linearized model ud         => analysis increment ua"
    endif

  elseif (xana%fields(index_ana)%fv3jedi_name == 'va') then

    !Special case: A-grid analysis, D-Grid model
    if (xmod%has_field('ud')) then
      !Already done above
      failed = .false.
      if (xana%f_comm%rank() == 0) write(*,"(A)") &
          "A2M MultiplyAdjoint: linearized model vd         => analysis increment va"
    endif

  elseif (xana%fields(index_ana)%fv3jedi_name == 'ps') then

    !Special case: ps in analysis, delp in model
    if (xmod%has_field('delp')) then
      call xana%get_field('ps',   xana_ps)
      call xmod%get_field('delp', xmod_delp)
      xana_ps = 0.0_kind_real
      do k = 1,geom%npz
        xana_ps(:,:,1) = xana_ps(:,:,1) + (geom%bk(k+1)-geom%bk(k))*xmod_delp(:,:,k)
      enddo
      failed = .false.
      if (xana%f_comm%rank() == 0) write(*,"(A)") &
          "A2M MultiplyAdjoint: linearized model delp       => analysis increment ps"
    endif

  endif

  if (failed) call abor1_ftn("fv3jedi_linvarcha_a2m_mod.multiplyadjoint: found no way of getting "//&
                             trim(xana%fields(index_ana)%fv3jedi_name)//" from the linearized model" )

enddo

end subroutine multiplyadjoint

! ------------------------------------------------------------------------------

subroutine multiplyinverse(self,geom,xmod,xana)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: xmod
type(fv3jedi_increment),     intent(inout) :: xana

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
        "A2M MultiplyInverse: linearized model "//xmod%fields(index_mod_found)%fv3jedi_name(1:10)&
        //" => analysis increment "//xana%fields(index_ana)%fv3jedi_name(1:10)

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
          "A2M MultiplyInverse: linearized model ud         => analysis increment ua"
    endif

  elseif (xana%fields(index_ana)%fv3jedi_name == 'va') then

    !Special case: A-grid analysis, D-Grid model
    if (xmod%has_field('ud')) then
      !Already done above
      failed = .false.
      if (xana%f_comm%rank() == 0) write(*,"(A)") &
          "A2M MultiplyInverse: linearized model vd         => analysis increment va"
    endif

  elseif (xana%fields(index_ana)%fv3jedi_name == 'ps') then

    !Special case: ps in analysis, delp in model
    if (xmod%has_field('delp')) then
      call xana%get_field('ps',   xana_ps)
      call xmod%get_field('delp', xmod_delp)
      xana_ps(:,:,1)   = sum(xmod_delp,3)
      failed = .false.
      if (xana%f_comm%rank() == 0) write(*,"(A)") &
          "A2M MultiplyInverse: linearized model delp       => analysis increment ps"
    endif

  endif

  if (failed) call abor1_ftn("fv3jedi_linvarcha_a2m_mod.multiplyinverse: found no way of getting "//&
                             trim(xana%fields(index_ana)%fv3jedi_name)//" from the linearized model" )

enddo

end subroutine multiplyinverse

! ------------------------------------------------------------------------------

subroutine multiplyinverseadjoint(self,geom,xana,xmod)

implicit none
type(fv3jedi_linvarcha_a2m), intent(inout) :: self
type(fv3jedi_geom),          intent(inout) :: geom
type(fv3jedi_increment),     intent(in)    :: xana
type(fv3jedi_increment),     intent(inout) :: xmod

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
    failed = .false.
    xmod%fields(index_mod)%array = xana%fields(index_ana_found)%array
    if (xmod%f_comm%rank() == 0) write(*,"(A, A10, A, A10)") &
        "A2M MultiplyInverseAdjoint: analysis increment "//xana%fields(index_ana_found)%fv3jedi_name(1:10)&
        //" => linearized model "//xmod%fields(index_mod)%fv3jedi_name(1:10)

  elseif (xmod%fields(index_mod)%fv3jedi_name == 'ud') then

    !Special case: A-grid analysis, D-Grid model
    if (xana%has_field('ua')) then
      call xana%get_field('ua', xana_ua)
      call xana%get_field('va', xana_va)
      call xmod%get_field('ud', xmod_ud)
      call xmod%get_field('vd', xmod_vd)
      call d_to_a_ad(geom, xmod_ud, xmod_vd, xana_ua, xana_va)
      xmod_ud(:,geom%jec+1,:) = 0.0_kind_real
      xmod_vd(geom%iec+1,:,:) = 0.0_kind_real
      failed = .false.
      if (xmod%f_comm%rank() == 0) write(*,"(A)") &
          "A2M MultiplyInverseAdjoint: analysis increment ua         => linearized model ud"
    endif

  elseif (xmod%fields(index_mod)%fv3jedi_name == 'vd') then

    !Special case: A-grid analysis, D-Grid model
    if (xana%has_field('ua')) then
      !Already done above
      failed = .false.
      if (xmod%f_comm%rank() == 0) write(*,"(A)") &
          "A2M MultiplyInverseAdjoint: analysis increment va         => linearized model vd"
    endif

  elseif (xmod%fields(index_mod)%fv3jedi_name == 'delp') then

    !Special case: ps in analysis, delp in model
    if (xana%has_field('ps')) then
      call xana%get_field('ps',   xana_ps)
      call xmod%get_field('delp', xmod_delp)
      do k = 1,geom%npz
        xmod_delp(:,:,k) = xana_ps(:,:,1)
      enddo
      failed = .false.
      if (xmod%f_comm%rank() == 0) write(*,"(A)") &
          "A2M MultiplyInverseAdjoint: analysis increment ps         => linearized model delp"
    endif

  endif

  if (failed) call abor1_ftn("fv3jedi_linvarcha_a2m_mod.multiplyinverseadjoint: found no way of getting "//&
                             trim(xmod%fields(index_mod)%fv3jedi_name)//" from the analysis increment" )

enddo

end subroutine multiplyinverseadjoint

! ------------------------------------------------------------------------------

end module fv3jedi_linvarcha_a2m_mod
