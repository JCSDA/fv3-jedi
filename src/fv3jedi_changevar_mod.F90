! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_changevar_mod

use fv3jedi_fields_mod, only: fv3jedi_field
use fv3jedi_geom_mod,   only: fv3jedi_geom
use iso_c_binding
use config_mod
use kinds

use moisture_vt_mod, only: esinit

implicit none

!> Fortran derived type to hold configuration data for the B mat variable change
type :: fv3jedi_changevar
 integer :: degsubs   = 100
 real(8) :: tmintbl   = 150.0, tmaxtbl = 333.0
 integer :: tablesize
 real(8), allocatable :: estblx(:)
 real(kind=kind_real), allocatable :: ttraj, qtraj, qsattraj
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

!Create lookup table for computing saturation specific humidity
self%tablesize = nint(self%tmaxtbl-self%tmintbl)*self%degsubs + 1
allocate(self%estblx(self%tablesize))
call esinit(self%tablesize,self%degsubs,self%tmintbl,self%tmaxtbl,self%estblx)

end subroutine fv3jedi_changevar_setup

! ------------------------------------------------------------------------------

subroutine fv3jedi_changevar_delete(self)

implicit none
type(fv3jedi_changevar), intent(inout) :: self

if (allocated(self%estblx)) deallocate(self%estblx)
if (allocated(self%ttraj)) deallocate(self%ttraj)
if (allocated(self%qtraj)) deallocate(self%qtraj)
if (allocated(self%qsattraj)) deallocate(self%qsattraj)

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
