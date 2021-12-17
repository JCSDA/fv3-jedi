! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! ------------------------------------------------------------------------------

module fv3jedi_linvarcha_nmcbal_interface_mod

use iso_c_binding

use fckit_configuration_module, only: fckit_configuration

use fv3jedi_linvarcha_nmcbal_mod
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_geom_interface_mod, only: fv3jedi_geom_registry
use fv3jedi_increment_interface_mod, only: fv3jedi_increment_registry
use fv3jedi_increment_mod, only: fv3jedi_increment

implicit none
private
public :: fv3jedi_linvarcha_nmcbal_registry

! ------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_linvarcha_nmcbal

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_linvarcha_nmcbal_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_nmcbal_create(c_key_self,c_key_geom,c_key_bg,c_key_fg,c_conf) &
           bind (c,name='fv3jedi_linvarcha_nmcbal_create_f90')

implicit none
integer(c_int),     intent(inout) :: c_key_self
integer(c_int),     intent(in)    :: c_key_geom
integer(c_int),     intent(in)    :: c_key_bg
integer(c_int),     intent(in)    :: c_key_fg
type(c_ptr), value, intent(in)    :: c_conf

type(fv3jedi_linvarcha_nmcbal), pointer :: self
type(fv3jedi_geom),             pointer :: geom
type(fv3jedi_state),            pointer :: bg
type(fv3jedi_state),            pointer :: fg
type(fckit_configuration) :: conf

call fv3jedi_linvarcha_nmcbal_registry%init()
call fv3jedi_linvarcha_nmcbal_registry%add(c_key_self)
call fv3jedi_linvarcha_nmcbal_registry%get(c_key_self, self)

call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_state_registry%get(c_key_bg,bg)
call fv3jedi_state_registry%get(c_key_fg,fg)

conf = fckit_configuration(c_conf)

call create(self,geom,bg,fg,conf)

end subroutine c_fv3jedi_linvarcha_nmcbal_create

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_nmcbal_delete(c_key_self) &
           bind (c,name='fv3jedi_linvarcha_nmcbal_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Change variable structure

type(fv3jedi_linvarcha_nmcbal), pointer :: self

call fv3jedi_linvarcha_nmcbal_registry%get(c_key_self,self)
call delete(self)
call fv3jedi_linvarcha_nmcbal_registry%remove(c_key_self)

end subroutine c_fv3jedi_linvarcha_nmcbal_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_nmcbal_multiply(c_key_self,c_key_geom,c_key_xuba,c_key_xbal) &
           bind (c,name='fv3jedi_linvarcha_nmcbal_multiply_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_xuba
integer(c_int), intent(in) :: c_key_xbal

type(fv3jedi_linvarcha_nmcbal), pointer :: self
type(fv3jedi_geom),             pointer :: geom
type(fv3jedi_increment),        pointer :: xuba
type(fv3jedi_increment),        pointer :: xbal

call fv3jedi_linvarcha_nmcbal_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_increment_registry%get(c_key_xuba,xuba)
call fv3jedi_increment_registry%get(c_key_xbal,xbal)

call multiply(self,geom,xuba,xbal)

end subroutine c_fv3jedi_linvarcha_nmcbal_multiply

! ----------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_nmcbal_multiplyadjoint(c_key_self,c_key_geom,c_key_xbal,c_key_xuba) &
           bind (c,name='fv3jedi_linvarcha_nmcbal_multiplyadjoint_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_xbal
integer(c_int), intent(in) :: c_key_xuba

type(fv3jedi_linvarcha_nmcbal), pointer :: self
type(fv3jedi_geom),             pointer :: geom
type(fv3jedi_increment),        pointer :: xbal
type(fv3jedi_increment),        pointer :: xuba

call fv3jedi_linvarcha_nmcbal_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_increment_registry%get(c_key_xbal,xbal)
call fv3jedi_increment_registry%get(c_key_xuba,xuba)

call multiplyadjoint(self,geom,xbal,xuba)

end subroutine c_fv3jedi_linvarcha_nmcbal_multiplyadjoint

! ----------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_nmcbal_multiplyinverse(c_key_self,c_key_geom,c_key_xbal,c_key_xuba) &
           bind (c,name='fv3jedi_linvarcha_nmcbal_multiplyinverse_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_xuba
integer(c_int), intent(in) :: c_key_xbal

type(fv3jedi_linvarcha_nmcbal), pointer :: self
type(fv3jedi_geom),             pointer :: geom
type(fv3jedi_increment),        pointer :: xbal
type(fv3jedi_increment),        pointer :: xuba

call fv3jedi_linvarcha_nmcbal_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_increment_registry%get(c_key_xbal,xbal)
call fv3jedi_increment_registry%get(c_key_xuba,xuba)

call multiplyinverse(self,geom,xbal,xuba)

end subroutine c_fv3jedi_linvarcha_nmcbal_multiplyinverse

! ----------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_nmcbal_multiplyinverseadjoint(c_key_self,c_key_geom,c_key_xuba,c_key_xbal) &
           bind (c,name='fv3jedi_linvarcha_nmcbal_multiplyinverseadjoint_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_xuba
integer(c_int), intent(in) :: c_key_xbal

type(fv3jedi_linvarcha_nmcbal), pointer :: self
type(fv3jedi_geom),             pointer :: geom
type(fv3jedi_increment),        pointer :: xbal
type(fv3jedi_increment),        pointer :: xuba

call fv3jedi_linvarcha_nmcbal_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_increment_registry%get(c_key_xuba,xuba)
call fv3jedi_increment_registry%get(c_key_xbal,xbal)

call multiplyinverseadjoint(self,geom,xuba,xbal)

end subroutine c_fv3jedi_linvarcha_nmcbal_multiplyinverseadjoint

! ----------------------------------------------------------------------------

end module fv3jedi_linvarcha_nmcbal_interface_mod
