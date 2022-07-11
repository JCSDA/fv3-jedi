! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! ------------------------------------------------------------------------------

module fv3jedi_linvarcha_a2m_interface_mod

use iso_c_binding
use fv3jedi_linvarcha_a2m_mod
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_geom_interface_mod, only: fv3jedi_geom_registry
use fv3jedi_increment_interface_mod, only: fv3jedi_increment_registry
use fv3jedi_increment_mod, only: fv3jedi_increment

implicit none
private
public :: fv3jedi_linvarcha_a2m_registry

! ------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_linvarcha_a2m

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_linvarcha_a2m_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_a2m_create(c_key_self,c_key_geom,c_key_bg,c_key_fg,c_conf) &
           bind (c,name='fv3jedi_linvarcha_a2m_create_f90')

implicit none
integer(c_int),     intent(inout) :: c_key_self
integer(c_int),     intent(in)    :: c_key_geom
integer(c_int),     intent(in)    :: c_key_bg
integer(c_int),     intent(in)    :: c_key_fg
type(c_ptr), value, intent(in)    :: c_conf

type(fv3jedi_linvarcha_a2m), pointer :: self
type(fv3jedi_geom),          pointer :: geom
type(fv3jedi_state),         pointer :: bg
type(fv3jedi_state),         pointer :: fg

call fv3jedi_linvarcha_a2m_registry%init()
call fv3jedi_linvarcha_a2m_registry%add(c_key_self)
call fv3jedi_linvarcha_a2m_registry%get(c_key_self, self)

call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_state_registry%get(c_key_bg,bg)
call fv3jedi_state_registry%get(c_key_fg,fg)

call create(self,geom,bg,fg,c_conf)

end subroutine c_fv3jedi_linvarcha_a2m_create

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_a2m_delete(c_key_self) &
           bind (c,name='fv3jedi_linvarcha_a2m_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Change variable structure

type(fv3jedi_linvarcha_a2m), pointer :: self

call fv3jedi_linvarcha_a2m_registry%get(c_key_self,self)
call delete(self)
call fv3jedi_linvarcha_a2m_registry%remove(c_key_self)

end subroutine c_fv3jedi_linvarcha_a2m_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_a2m_multiply(c_key_self,c_key_geom,c_key_xana,c_key_xmod) &
           bind (c,name='fv3jedi_linvarcha_a2m_multiply_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_xana
integer(c_int), intent(in) :: c_key_xmod

type(fv3jedi_linvarcha_a2m), pointer :: self
type(fv3jedi_geom),          pointer :: geom
type(fv3jedi_increment),     pointer :: xana
type(fv3jedi_increment),     pointer :: xmod

call fv3jedi_linvarcha_a2m_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_increment_registry%get(c_key_xana,xana)
call fv3jedi_increment_registry%get(c_key_xmod,xmod)

call multiply(self,geom,xana,xmod)

end subroutine c_fv3jedi_linvarcha_a2m_multiply

! ----------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_a2m_multiplyadjoint(c_key_self,c_key_geom,c_key_xmod,c_key_xana) &
           bind (c,name='fv3jedi_linvarcha_a2m_multiplyadjoint_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_xmod
integer(c_int), intent(in) :: c_key_xana

type(fv3jedi_linvarcha_a2m), pointer :: self
type(fv3jedi_geom),          pointer :: geom
type(fv3jedi_increment),     pointer :: xmod
type(fv3jedi_increment),     pointer :: xana

call fv3jedi_linvarcha_a2m_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_increment_registry%get(c_key_xmod,xmod)
call fv3jedi_increment_registry%get(c_key_xana,xana)

call multiplyadjoint(self,geom,xmod,xana)

end subroutine c_fv3jedi_linvarcha_a2m_multiplyadjoint

! ----------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_a2m_multiplyinverse(c_key_self,c_key_geom,c_key_xmod,c_key_xana) &
           bind (c,name='fv3jedi_linvarcha_a2m_multiplyinverse_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_xmod
integer(c_int), intent(in) :: c_key_xana

type(fv3jedi_linvarcha_a2m), pointer :: self
type(fv3jedi_geom),          pointer :: geom
type(fv3jedi_increment),     pointer :: xmod
type(fv3jedi_increment),     pointer :: xana

call fv3jedi_linvarcha_a2m_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_increment_registry%get(c_key_xmod,xmod)
call fv3jedi_increment_registry%get(c_key_xana,xana)

call multiplyinverse(self,geom,xmod,xana)

end subroutine c_fv3jedi_linvarcha_a2m_multiplyinverse

! ----------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_a2m_multiplyinverseadjoint(c_key_self,c_key_geom,c_key_xana,c_key_xmod) &
           bind (c,name='fv3jedi_linvarcha_a2m_multiplyinverseadjoint_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_xana
integer(c_int), intent(in) :: c_key_xmod

type(fv3jedi_linvarcha_a2m), pointer :: self
type(fv3jedi_geom),          pointer :: geom
type(fv3jedi_increment),     pointer :: xana
type(fv3jedi_increment),     pointer :: xmod

call fv3jedi_linvarcha_a2m_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_increment_registry%get(c_key_xana,xana)
call fv3jedi_increment_registry%get(c_key_xmod,xmod)

call multiplyinverseadjoint(self,geom,xana,xmod)

end subroutine c_fv3jedi_linvarcha_a2m_multiplyinverseadjoint

! ----------------------------------------------------------------------------

end module fv3jedi_linvarcha_a2m_interface_mod
