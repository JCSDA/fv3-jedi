! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! ------------------------------------------------------------------------------

module fv3jedi_vc_geosrst2bkg_interface_mod

use iso_c_binding
use datetime_mod

use fckit_configuration_module, only: fckit_configuration

use fv3jedi_vc_geosrst2bkg_mod
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_geom_interface_mod, only: fv3jedi_geom_registry

implicit none
private
public :: fv3jedi_vc_geosrst2bkg_registry

! ------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_vc_geosrst2bkg

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_vc_geosrst2bkg_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_vc_geosrst2bkg_create(c_key_self,c_key_geom,c_conf) &
           bind (c,name='fv3jedi_vc_geosrst2bkg_create_f90')

implicit none
integer(c_int),     intent(inout) :: c_key_self
integer(c_int),     intent(in)    :: c_key_geom
type(c_ptr), value, intent(in)    :: c_conf

type(fv3jedi_vc_geosrst2bkg), pointer :: self
type(fv3jedi_geom),           pointer :: geom
type(fckit_configuration)             :: conf

call fv3jedi_vc_geosrst2bkg_registry%init()
call fv3jedi_vc_geosrst2bkg_registry%add(c_key_self)
call fv3jedi_vc_geosrst2bkg_registry%get(c_key_self, self)

call fv3jedi_geom_registry%get(c_key_geom,geom)

conf = fckit_configuration(c_conf)

call create(self,geom,conf)

end subroutine c_fv3jedi_vc_geosrst2bkg_create

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_vc_geosrst2bkg_delete(c_key_self) &
           bind (c,name='fv3jedi_vc_geosrst2bkg_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Change variable structure

type(fv3jedi_vc_geosrst2bkg), pointer :: self

call fv3jedi_vc_geosrst2bkg_registry%get(c_key_self,self)
call delete(self)
call fv3jedi_vc_geosrst2bkg_registry%remove(c_key_self)

end subroutine c_fv3jedi_vc_geosrst2bkg_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_vc_geosrst2bkg_changevar(c_key_self,c_key_geom,c_key_xd,c_key_xa) &
           bind (c,name='fv3jedi_vc_geosrst2bkg_changevar_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_xd
integer(c_int), intent(in) :: c_key_xa

type(fv3jedi_vc_geosrst2bkg), pointer :: self
type(fv3jedi_geom),       pointer :: geom
type(fv3jedi_state),      pointer :: xr
type(fv3jedi_state),      pointer :: xb

call fv3jedi_vc_geosrst2bkg_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_state_registry%get(c_key_xd,xr)
call fv3jedi_state_registry%get(c_key_xa,xb)

call changevar(self,geom,xr,xb)

end subroutine c_fv3jedi_vc_geosrst2bkg_changevar

! ----------------------------------------------------------------------------

subroutine c_fv3jedi_vc_geosrst2bkg_changevarinverse(c_key_self,c_key_geom,c_key_xa,c_key_xd) &
           bind (c,name='fv3jedi_vc_geosrst2bkg_changevarinverse_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_xa
integer(c_int), intent(in) :: c_key_xd

type(fv3jedi_vc_geosrst2bkg), pointer :: self
type(fv3jedi_geom),       pointer :: geom
type(fv3jedi_state),      pointer :: xr
type(fv3jedi_state),      pointer :: xb

call fv3jedi_vc_geosrst2bkg_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_state_registry%get(c_key_xa,xb)
call fv3jedi_state_registry%get(c_key_xd,xr)

call changevarinverse(self,geom,xb,xr)

end subroutine c_fv3jedi_vc_geosrst2bkg_changevarinverse

! ----------------------------------------------------------------------------

end module fv3jedi_vc_geosrst2bkg_interface_mod
