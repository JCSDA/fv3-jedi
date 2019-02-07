! (C) Copyright 2018 UCAR
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
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_linvarcha_a2m_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_a2m_setup(c_key_self, c_key_geom, c_key_state_bg, c_key_state_fg, &
           c_conf) bind (c,name='fv3jedi_linvarcha_a2m_setup_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self     !< Change variable structure
integer(c_int), intent(in)    :: c_key_state_bg !< Background key
integer(c_int), intent(in)    :: c_key_state_fg !< First guess key
integer(c_int), intent(in)    :: c_key_geom     !< Geom key
type(c_ptr),    intent(in)    :: c_conf         !< Configuration

type(fv3jedi_linvarcha_a2m), pointer :: self
type(fv3jedi_state), pointer :: bg
type(fv3jedi_state), pointer :: fg
type(fv3jedi_geom), pointer :: geom

call fv3jedi_linvarcha_a2m_registry%init()
call fv3jedi_linvarcha_a2m_registry%add(c_key_self)
call fv3jedi_linvarcha_a2m_registry%get(c_key_self, self)

call fv3jedi_state_registry%get(c_key_state_bg,bg)
call fv3jedi_state_registry%get(c_key_state_fg,fg)

call fv3jedi_geom_registry%get(c_key_geom,geom)

call create(self, bg, fg, geom, c_conf)

end subroutine c_fv3jedi_linvarcha_a2m_setup

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

subroutine c_fv3jedi_linvarcha_a2m_multiply(c_key_self, c_key_geom, c_key_in, c_key_out) &
           bind (c,name='fv3jedi_linvarcha_a2m_multiply_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom     !< Geom key
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_linvarcha_a2m), pointer :: self
type(fv3jedi_increment), pointer :: xin
type(fv3jedi_increment), pointer :: xout
type(fv3jedi_geom), pointer :: geom

call fv3jedi_linvarcha_a2m_registry%get(c_key_self,self)
call fv3jedi_increment_registry%get(c_key_in,xin)
call fv3jedi_increment_registry%get(c_key_out,xout)
call fv3jedi_geom_registry%get(c_key_geom,geom)

call multiply(self,geom,xin,xout)

end subroutine c_fv3jedi_linvarcha_a2m_multiply

! ----------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_a2m_multiplyadjoint(c_key_self, c_key_geom, c_key_in, &
           c_key_out) bind (c,name='fv3jedi_linvarcha_a2m_multiplyadjoint_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom     !< Geom key
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_linvarcha_a2m), pointer :: self
type(fv3jedi_increment), pointer :: xin
type(fv3jedi_increment), pointer :: xout
type(fv3jedi_geom), pointer :: geom

call fv3jedi_linvarcha_a2m_registry%get(c_key_self,self)
call fv3jedi_increment_registry%get(c_key_in,xin)
call fv3jedi_increment_registry%get(c_key_out,xout)
call fv3jedi_geom_registry%get(c_key_geom,geom)

call multiplyadjoint(self,geom,xin,xout)

end subroutine c_fv3jedi_linvarcha_a2m_multiplyadjoint

! ----------------------------------------------------------------------------

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_a2m_multiplyinverse(c_key_self, c_key_geom, c_key_in, c_key_out) &
           bind (c,name='fv3jedi_linvarcha_a2m_multiplyinverse_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom     !< Geom key
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_linvarcha_a2m), pointer :: self
type(fv3jedi_increment), pointer :: xin
type(fv3jedi_increment), pointer :: xout
type(fv3jedi_geom), pointer :: geom

call fv3jedi_linvarcha_a2m_registry%get(c_key_self,self)
call fv3jedi_increment_registry%get(c_key_in,xin)
call fv3jedi_increment_registry%get(c_key_out,xout)
call fv3jedi_geom_registry%get(c_key_geom,geom)

call multiplyinverse(self,geom,xin,xout)

end subroutine c_fv3jedi_linvarcha_a2m_multiplyinverse

! ----------------------------------------------------------------------------

subroutine c_fv3jedi_linvarcha_a2m_multiplyinverseadjoint(c_key_self, c_key_geom, c_key_in, &
      c_key_out) bind (c,name='fv3jedi_linvarcha_a2m_multiplyinverseadjoint_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom     !< Geom key
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_linvarcha_a2m), pointer :: self
type(fv3jedi_increment), pointer :: xin
type(fv3jedi_increment), pointer :: xout
type(fv3jedi_geom), pointer :: geom

call fv3jedi_linvarcha_a2m_registry%get(c_key_self,self)
call fv3jedi_increment_registry%get(c_key_in,xin)
call fv3jedi_increment_registry%get(c_key_out,xout)
call fv3jedi_geom_registry%get(c_key_geom,geom)

call multiplyinverseadjoint(self,geom,xin,xout)

end subroutine c_fv3jedi_linvarcha_a2m_multiplyinverseadjoint

! ----------------------------------------------------------------------------

end module fv3jedi_linvarcha_a2m_interface_mod
