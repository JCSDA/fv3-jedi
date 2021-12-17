! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! --------------------------------------------------------------------------------------------------

module fv3jedi_vc_coldstartwinds_interface_mod

use iso_c_binding
use datetime_mod

use fckit_configuration_module, only: fckit_configuration

use fv3jedi_vc_coldstartwinds_mod
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_geom_interface_mod, only: fv3jedi_geom_registry

implicit none
private
public :: fv3jedi_vc_coldstartwinds_registry

! --------------------------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_vc_coldstartwinds

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_vc_coldstartwinds_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_vc_coldstartwinds_create(c_key_self,c_key_geom,c_conf) &
           bind (c,name='fv3jedi_vc_coldstartwinds_create_f90')

integer(c_int),     intent(inout) :: c_key_self
integer(c_int),     intent(in)    :: c_key_geom
type(c_ptr), value, intent(in)    :: c_conf

type(fv3jedi_vc_coldstartwinds), pointer :: self
type(fv3jedi_geom),              pointer :: geom
type(fckit_configuration)                :: conf

! LinkedList
call fv3jedi_vc_coldstartwinds_registry%init()
call fv3jedi_vc_coldstartwinds_registry%add(c_key_self)
call fv3jedi_vc_coldstartwinds_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_geom,geom)

conf = fckit_configuration(c_conf)

! Implementation
call self%create(geom, conf)

end subroutine c_fv3jedi_vc_coldstartwinds_create

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_vc_coldstartwinds_delete(c_key_self) &
           bind (c,name='fv3jedi_vc_coldstartwinds_delete_f90')

integer(c_int), intent(inout) :: c_key_self  !< Change variable structure

type(fv3jedi_vc_coldstartwinds), pointer :: self

! LinkedList
call fv3jedi_vc_coldstartwinds_registry%get(c_key_self,self)

! Implementation
call self%delete()

! LinkedList
call fv3jedi_vc_coldstartwinds_registry%remove(c_key_self)

end subroutine c_fv3jedi_vc_coldstartwinds_delete

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_vc_coldstartwinds_changevar(c_key_self, c_key_xin, c_key_xout) &
           bind (c,name='fv3jedi_vc_coldstartwinds_changevar_f90')

integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_xin
integer(c_int), intent(in) :: c_key_xout

type(fv3jedi_vc_coldstartwinds), pointer :: self
type(fv3jedi_state),             pointer :: xin
type(fv3jedi_state),             pointer :: xout

! LinkedList
call fv3jedi_vc_coldstartwinds_registry%get(c_key_self,self)
call fv3jedi_state_registry%get(c_key_xin,xin)
call fv3jedi_state_registry%get(c_key_xout,xout)

! Implementation
call self%changevar(xin, xout)

end subroutine c_fv3jedi_vc_coldstartwinds_changevar

! --------------------------------------------------------------------------------------------------

end module fv3jedi_vc_coldstartwinds_interface_mod
