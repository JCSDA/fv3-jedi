! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_changevar_setup(c_key_self, c_conf) &
           bind (c,name='fv3jedi_changevar_setup_f90')

use iso_c_binding
use fv3jedi_changevar_mod

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Change variable structure
type(c_ptr),    intent(in)    :: c_conf      !< Configuration

type(fv3jedi_changevar), pointer :: self

call fv3jedi_changevar_registry%init()
call fv3jedi_changevar_registry%add(c_key_self)
call fv3jedi_changevar_registry%get(c_key_self, self)

call fv3jedi_changevar_setup(self, c_conf)

end subroutine c_fv3jedi_changevar_setup

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_changevar_delete(c_key_self) &
           bind (c,name='fv3jedi_changevar_delete_f90')

use iso_c_binding
use fv3jedi_changevar_mod

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Change variable structure

type(fv3jedi_changevar), pointer :: self

call fv3jedi_changevar_registry%get(c_key_self,self)
call fv3jedi_changevar_delete(self)
call fv3jedi_changevar_registry%remove(c_key_self)

end subroutine c_fv3jedi_changevar_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_changevar_linearize(c_key_self,c_key_geom,c_key_state) &
           bind (c,name='fv3jedi_changevar_linearize_f90')

use iso_c_binding
use fv3jedi_changevar_mod
use fv3jedi_geom_mod
use fv3jedi_fields_mod

implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(inout) :: c_key_geom
integer(c_int), intent(inout) :: c_key_state

type(fv3jedi_changevar), pointer :: self
type(fv3jedi_geom), pointer :: geom
type(fv3jedi_field), pointer :: state

call fv3jedi_changevar_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_field_registry%get(c_key_state,state)

call fv3jedi_changevar_linearize(self,geom,state)

end subroutine c_fv3jedi_changevar_linearize

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_changevar_transform(c_key_self, c_key_in, c_key_out) &
           bind (c,name='fv3jedi_changevar_transform_f90')

use iso_c_binding
use fv3jedi_changevar_mod
use fv3jedi_fields_mod, only: fv3jedi_field_registry, fv3jedi_field
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_changevar), pointer :: self
type(fv3jedi_field), pointer :: xin
type(fv3jedi_field), pointer :: xout

call fv3jedi_changevar_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_in,xin)
call fv3jedi_field_registry%get(c_key_out,xout)

call fv3jedi_changevar_transform(self,xin,xout)

end subroutine c_fv3jedi_changevar_transform

! ----------------------------------------------------------------------------

subroutine c_fv3jedi_changevar_transformadjoint(c_key_self, c_key_in, &
           c_key_out) bind (c,name='fv3jedi_changevar_transformadjoint_f90')

use iso_c_binding
use fv3jedi_changevar_mod
use fv3jedi_fields_mod, only: fv3jedi_field_registry, fv3jedi_field
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_changevar), pointer :: self
type(fv3jedi_field), pointer :: xin
type(fv3jedi_field), pointer :: xout

call fv3jedi_changevar_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_in,xin)
call fv3jedi_field_registry%get(c_key_out,xout)

call fv3jedi_changevar_transformadjoint(self,xin,xout)

end subroutine c_fv3jedi_changevar_transformadjoint

! ----------------------------------------------------------------------------

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_changevar_transforminverse(c_key_self, c_key_in, c_key_out) &
           bind (c,name='fv3jedi_changevar_transforminverse_f90')

use iso_c_binding
use fv3jedi_changevar_mod
use fv3jedi_fields_mod, only: fv3jedi_field_registry, fv3jedi_field
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_changevar), pointer :: self
type(fv3jedi_field), pointer :: xin
type(fv3jedi_field), pointer :: xout

call fv3jedi_changevar_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_in,xin)
call fv3jedi_field_registry%get(c_key_out,xout)

call fv3jedi_changevar_transforminverse(self,xin,xout)

end subroutine c_fv3jedi_changevar_transforminverse

! ----------------------------------------------------------------------------

subroutine c_fv3jedi_changevar_transforminverseadjoint(c_key_self, c_key_in, &
      c_key_out) bind (c,name='fv3jedi_changevar_transforminverseadjoint_f90')

use iso_c_binding
use fv3jedi_changevar_mod
use fv3jedi_fields_mod, only: fv3jedi_field_registry, fv3jedi_field
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_changevar), pointer :: self
type(fv3jedi_field), pointer :: xin
type(fv3jedi_field), pointer :: xout

call fv3jedi_changevar_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_in,xin)
call fv3jedi_field_registry%get(c_key_out,xout)

call fv3jedi_changevar_transforminverseadjoint(self,xin,xout)

end subroutine c_fv3jedi_changevar_transforminverseadjoint

! ----------------------------------------------------------------------------
