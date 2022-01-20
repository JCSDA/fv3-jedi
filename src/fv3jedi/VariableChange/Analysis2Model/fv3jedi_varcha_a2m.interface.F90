! (C) Copyright 2018-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! --------------------------------------------------------------------------------------------------

module fv3jedi_varcha_a2m_interface_mod

! iso
use iso_c_binding

! fckit
use fckit_configuration_module, only: fckit_configuration

! oops
use datetime_mod

! fv3-jedi
use fv3jedi_varcha_a2m_mod,      only: fv3jedi_varcha_a2m
use fv3jedi_state_mod,           only: fv3jedi_state
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_geom_interface_mod,  only: fv3jedi_geom_registry

implicit none
private
public :: fv3jedi_varcha_a2m_registry

! --------------------------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_varcha_a2m

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_varcha_a2m_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_varcha_a2m_create(c_key_self,c_key_geom,c_conf) &
           bind (c,name='fv3jedi_varcha_a2m_create_f90')

implicit none
integer(c_int),     intent(inout) :: c_key_self
integer(c_int),     intent(in)    :: c_key_geom
type(c_ptr), value, intent(in)    :: c_conf

type(fv3jedi_varcha_a2m), pointer :: self
type(fv3jedi_geom),       pointer :: geom
type(fckit_configuration)         :: conf

call fv3jedi_varcha_a2m_registry%init()
call fv3jedi_varcha_a2m_registry%add(c_key_self)
call fv3jedi_varcha_a2m_registry%get(c_key_self, self)

call fv3jedi_geom_registry%get(c_key_geom,geom)

! Fortran configuration
! ---------------------
conf = fckit_configuration(c_conf)

call self%create(geom, conf)

end subroutine c_fv3jedi_varcha_a2m_create

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_varcha_a2m_delete(c_key_self) &
           bind (c,name='fv3jedi_varcha_a2m_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Change variable structure

type(fv3jedi_varcha_a2m), pointer :: self

call fv3jedi_varcha_a2m_registry%get(c_key_self,self)
call self%delete()
call fv3jedi_varcha_a2m_registry%remove(c_key_self)

end subroutine c_fv3jedi_varcha_a2m_delete

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_varcha_a2m_changevar(c_key_self,c_key_geom,c_key_xana,c_key_xmod) &
           bind (c,name='fv3jedi_varcha_a2m_changevar_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_xana
integer(c_int), intent(in) :: c_key_xmod

type(fv3jedi_varcha_a2m), pointer :: self
type(fv3jedi_geom),       pointer :: geom
type(fv3jedi_state),      pointer :: xana
type(fv3jedi_state),      pointer :: xmod

call fv3jedi_varcha_a2m_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_state_registry%get(c_key_xana,xana)
call fv3jedi_state_registry%get(c_key_xmod,xmod)

call self%changevar(geom,xana,xmod)

end subroutine c_fv3jedi_varcha_a2m_changevar

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_varcha_a2m_changevarinverse(c_key_self,c_key_geom,c_key_xmod,c_key_xana) &
           bind (c,name='fv3jedi_varcha_a2m_changevarinverse_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_xmod
integer(c_int), intent(in) :: c_key_xana

type(fv3jedi_varcha_a2m), pointer :: self
type(fv3jedi_geom),       pointer :: geom
type(fv3jedi_state),      pointer :: xana
type(fv3jedi_state),      pointer :: xmod

call fv3jedi_varcha_a2m_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_state_registry%get(c_key_xmod,xmod)
call fv3jedi_state_registry%get(c_key_xana,xana)

call self%changevarinverse(geom,xmod,xana)

end subroutine c_fv3jedi_varcha_a2m_changevarinverse

! --------------------------------------------------------------------------------------------------

end module fv3jedi_varcha_a2m_interface_mod
