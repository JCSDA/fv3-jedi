! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! --------------------------------------------------------------------------------------------------

module fv3jedi_vc_model2geovals_interface_mod

use iso_c_binding

use datetime_mod

use fckit_configuration_module, only: fckit_configuration

use fv3jedi_geom_mod,             only: fv3jedi_geom
use fv3jedi_geom_interface_mod,   only: fv3jedi_geom_registry
use fv3jedi_state_mod,            only: fv3jedi_state
use fv3jedi_state_interface_mod,  only: fv3jedi_state_registry
use fv3jedi_vc_model2geovals_mod, only: fv3jedi_vc_model2geovals

implicit none

private
public :: fv3jedi_vc_model2geovals_registry

! --------------------------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_vc_model2geovals

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_vc_model2geovals_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_vc_model2geovals_create(c_key_self, c_key_geom, c_conf) &
           bind (c, name='fv3jedi_vc_model2geovals_create_f90')

implicit none
integer(c_int),     intent(inout) :: c_key_self
integer(c_int),     intent(in)    :: c_key_geom
type(c_ptr), value, intent(in)    :: c_conf

type(fv3jedi_vc_model2geovals), pointer :: self
type(fv3jedi_geom),             pointer :: geom
type(fckit_configuration)               :: conf

! Linked list
! -----------
call fv3jedi_vc_model2geovals_registry%init()
call fv3jedi_vc_model2geovals_registry%add(c_key_self)
call fv3jedi_vc_model2geovals_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_geom,geom)

! APIs
! ----
conf = fckit_configuration(c_conf)

! Implementation
! --------------
call self%create(geom, conf)

end subroutine c_fv3jedi_vc_model2geovals_create

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_vc_model2geovals_delete(c_key_self) &
           bind (c, name='fv3jedi_vc_model2geovals_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Change variable structure

type(fv3jedi_vc_model2geovals), pointer :: self

! Linked list
! -----------
call fv3jedi_vc_model2geovals_registry%get(c_key_self,self)

! Implementation
! --------------
call self%delete()

! Linked list
! -----------
call fv3jedi_vc_model2geovals_registry%remove(c_key_self)

end subroutine c_fv3jedi_vc_model2geovals_delete

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_vc_model2geovals_changevar(c_key_self, c_key_geom, c_key_xm, c_key_xg) &
           bind (c, name='fv3jedi_vc_model2geovals_changevar_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_xm
integer(c_int), intent(in) :: c_key_xg

type(fv3jedi_vc_model2geovals), pointer :: self
type(fv3jedi_geom),             pointer :: geom
type(fv3jedi_state),            pointer :: xm
type(fv3jedi_state),            pointer :: xg

! Linked list
! -----------
call fv3jedi_vc_model2geovals_registry%get(c_key_self,self)
call fv3jedi_state_registry%get(c_key_xm,xm)
call fv3jedi_state_registry%get(c_key_xg,xg)
call fv3jedi_geom_registry%get(c_key_geom,geom)

! Implementation
! --------------
call self%changevar(geom, xm, xg)

end subroutine c_fv3jedi_vc_model2geovals_changevar

! --------------------------------------------------------------------------------------------------

end module fv3jedi_vc_model2geovals_interface_mod
