! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_changevar_mod

use fv3jedi_fields_mod, only: fv3jedi_field
use fv3jedi_geom_mod,   only: fv3jedi_geom
use iso_c_binding
use config_mod

implicit none

!> Fortran derived type to hold configuration data for the B mat variable change
type :: fv3jedi_changevar
  integer :: nothing_yet
end type fv3jedi_changevar

#define LISTED_TYPE fv3jedi_changevar

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_changevar_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "linkedList_c.f"
! ------------------------------------------------------------------------------

subroutine fv3jedi_changevar_setup(self, c_conf)

implicit none
type(fv3jedi_changevar), intent(inout) :: self    !< Change variable structure
type(c_ptr),             intent(in)    :: c_conf  !< Configuration

end subroutine fv3jedi_changevar_setup

! ------------------------------------------------------------------------------

subroutine fv3jedi_changevar_delete(self)

implicit none
type(fv3jedi_changevar), intent(inout) :: self

end subroutine fv3jedi_changevar_delete

! ------------------------------------------------------------------------------

subroutine fv3jedi_changevar_linearize(self,geom,state)

implicit none
type(fv3jedi_changevar), intent(inout) :: self
type(fv3jedi_geom),  intent(inout) :: geom
type(fv3jedi_field), intent(inout) :: state

end subroutine fv3jedi_changevar_linearize

! ------------------------------------------------------------------------------

subroutine fv3jedi_changevar_transform(self,xinc,xctr)

implicit none
type(fv3jedi_changevar), intent(inout) :: self
type(fv3jedi_field), intent(inout) :: xinc
type(fv3jedi_field), intent(inout) :: xctr

end subroutine fv3jedi_changevar_transform

! ------------------------------------------------------------------------------

subroutine fv3jedi_changevar_transformadjoint(self,xinc,xctr)

implicit none
type(fv3jedi_changevar), intent(inout) :: self
type(fv3jedi_field), intent(inout) :: xinc
type(fv3jedi_field), intent(inout) :: xctr

end subroutine fv3jedi_changevar_transformadjoint

! ------------------------------------------------------------------------------

subroutine fv3jedi_changevar_transforminverse(self,xinc,xctr)

implicit none
type(fv3jedi_changevar), intent(inout) :: self
type(fv3jedi_field), intent(inout) :: xinc
type(fv3jedi_field), intent(inout) :: xctr

end subroutine fv3jedi_changevar_transforminverse

! ------------------------------------------------------------------------------

subroutine fv3jedi_changevar_transforminverseadjoint(self,xinc,xctr)

implicit none
type(fv3jedi_changevar), intent(inout) :: self
type(fv3jedi_field), intent(inout) :: xinc
type(fv3jedi_field), intent(inout) :: xctr

end subroutine fv3jedi_changevar_transforminverseadjoint

! ------------------------------------------------------------------------------

end module fv3jedi_changevar_mod
