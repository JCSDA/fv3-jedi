! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_lvc_model2geovals_interface_mod

use iso_c_binding

use datetime_mod

use fckit_configuration_module, only: fckit_configuration

use fv3jedi_geom_mod,                only: fv3jedi_geom
use fv3jedi_geom_interface_mod,      only: fv3jedi_geom_registry
use fv3jedi_increment_interface_mod, only: fv3jedi_increment_registry
use fv3jedi_increment_mod,           only: fv3jedi_increment
use fv3jedi_lvc_model2geovals_mod,   only: fv3jedi_lvc_model2geovals
use fv3jedi_state_interface_mod,     only: fv3jedi_state_registry
use fv3jedi_state_mod,               only: fv3jedi_state

implicit none

private
public :: fv3jedi_lvc_model2geovals_registry

! --------------------------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_lvc_model2geovals

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_lvc_model2geovals_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_lvc_model2geovals_create(c_key_self, c_key_geom, c_key_bg, c_key_fg, c_conf) &
           bind (c,name='fv3jedi_lvc_model2geovals_create_f90')

integer(c_int),     intent(inout) :: c_key_self
integer(c_int),     intent(in)    :: c_key_geom
integer(c_int),     intent(in)    :: c_key_bg
integer(c_int),     intent(in)    :: c_key_fg
type(c_ptr), value, intent(in)    :: c_conf

type(fv3jedi_lvc_model2geovals), pointer :: self
type(fv3jedi_geom),              pointer :: geom
type(fv3jedi_state),             pointer :: bg
type(fv3jedi_state),             pointer :: fg
type(fckit_configuration)                :: conf

! Linked list
! -----------
call fv3jedi_lvc_model2geovals_registry%init()
call fv3jedi_lvc_model2geovals_registry%add(c_key_self)
call fv3jedi_lvc_model2geovals_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_state_registry%get(c_key_bg, bg)
call fv3jedi_state_registry%get(c_key_fg, fg)

! APIs
! ----
conf = fckit_configuration(c_conf)

! Implementation
! --------------
call self%create(geom, bg, fg, conf)

end subroutine c_fv3jedi_lvc_model2geovals_create

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_lvc_model2geovals_delete(c_key_self) &
           bind (c,name='fv3jedi_lvc_model2geovals_delete_f90')

integer(c_int), intent(inout) :: c_key_self

type(fv3jedi_lvc_model2geovals), pointer :: self

! Linked list
! -----------
call fv3jedi_lvc_model2geovals_registry%get(c_key_self, self)

! Implementation
! --------------
call self%delete()

! Linked list
! -----------
call fv3jedi_lvc_model2geovals_registry%remove(c_key_self)

end subroutine c_fv3jedi_lvc_model2geovals_delete

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_lvc_model2geovals_multiply(c_key_self, c_key_geom, c_key_dxm, c_key_dxg) &
           bind (c,name='fv3jedi_lvc_model2geovals_multiply_f90')

integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_dxm
integer(c_int), intent(in) :: c_key_dxg

type(fv3jedi_lvc_model2geovals), pointer :: self
type(fv3jedi_geom),              pointer :: geom
type(fv3jedi_increment),         pointer :: dxm
type(fv3jedi_increment),         pointer :: dxg

! Linked list
! -----------
call fv3jedi_lvc_model2geovals_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_increment_registry%get(c_key_dxm,dxm)
call fv3jedi_increment_registry%get(c_key_dxg,dxg)

! Implementation
! --------------
call self%multiply(geom, dxm, dxg)

end subroutine c_fv3jedi_lvc_model2geovals_multiply

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_lvc_model2geovals_multiplyadjoint(c_key_self, c_key_geom, c_key_dxg, &
                                                       c_key_dxm) &
           bind (c,name='fv3jedi_lvc_model2geovals_multiplyadjoint_f90')

integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_dxg
integer(c_int), intent(in) :: c_key_dxm

type(fv3jedi_lvc_model2geovals), pointer :: self
type(fv3jedi_geom),              pointer :: geom
type(fv3jedi_increment),         pointer :: dxg
type(fv3jedi_increment),         pointer :: dxm

! Linked list
! -----------
call fv3jedi_lvc_model2geovals_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_increment_registry%get(c_key_dxg, dxg)
call fv3jedi_increment_registry%get(c_key_dxm, dxm)

! Implementation
! --------------
call self%multiplyadjoint(geom, dxg, dxm)

end subroutine c_fv3jedi_lvc_model2geovals_multiplyadjoint

! --------------------------------------------------------------------------------------------------

end module fv3jedi_lvc_model2geovals_interface_mod
