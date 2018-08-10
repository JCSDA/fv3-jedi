! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_varchange_setup(c_key_self, c_conf) &
           bind (c,name='fv3jedi_varchange_setup_f90')

use iso_c_binding
use fv3jedi_varchange_mod

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Change variable structure
type(c_ptr),    intent(in)    :: c_conf      !< Configuration

type(fv3jedi_varchange), pointer :: self

call fv3jedi_varchange_registry%init()
call fv3jedi_varchange_registry%add(c_key_self)
call fv3jedi_varchange_registry%get(c_key_self, self)

call fv3jedi_varchange_setup(self, c_conf)

end subroutine c_fv3jedi_varchange_setup

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_varchange_delete(c_key_self) &
           bind (c,name='fv3jedi_varchange_delete_f90')

use iso_c_binding
use fv3jedi_varchange_mod

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Change variable structure

type(fv3jedi_varchange), pointer :: self

call fv3jedi_varchange_registry%get(c_key_self,self)
call fv3jedi_varchange_delete(self)
call fv3jedi_varchange_registry%remove(c_key_self)

end subroutine c_fv3jedi_varchange_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_varchange_linearize(c_key_self,c_key_geom,c_key_state) &
           bind (c,name='fv3jedi_varchange_linearize_f90')

use iso_c_binding
use fv3jedi_varchange_mod
use fv3jedi_geom_mod
use fv3jedi_fields_mod

implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(inout) :: c_key_geom
integer(c_int), intent(inout) :: c_key_state

type(fv3jedi_varchange), pointer :: self
type(fv3jedi_geom), pointer :: geom
type(fv3jedi_field), pointer :: state

call fv3jedi_varchange_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_field_registry%get(c_key_state,state)

call fv3jedi_varchange_linearize(self,geom,state)

end subroutine c_fv3jedi_varchange_linearize

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_varchange_multiply(c_key_self, c_key_in, c_key_out) &
           bind (c,name='fv3jedi_varchange_multiply_f90')

use iso_c_binding
use fv3jedi_varchange_mod
use fv3jedi_fields_mod, only: fv3jedi_field_registry, fv3jedi_field
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_varchange), pointer :: self
type(fv3jedi_field), pointer :: xin
type(fv3jedi_field), pointer :: xout

call fv3jedi_varchange_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_in,xin)
call fv3jedi_field_registry%get(c_key_out,xout)

call fv3jedi_varchange_multiply(self,xin,xout)

end subroutine c_fv3jedi_varchange_multiply

! ----------------------------------------------------------------------------

subroutine c_fv3jedi_varchange_multiplyadjoint(c_key_self, c_key_in, &
           c_key_out) bind (c,name='fv3jedi_varchange_multiplyadjoint_f90')

use iso_c_binding
use fv3jedi_varchange_mod
use fv3jedi_fields_mod, only: fv3jedi_field_registry, fv3jedi_field
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_varchange), pointer :: self
type(fv3jedi_field), pointer :: xin
type(fv3jedi_field), pointer :: xout

call fv3jedi_varchange_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_in,xin)
call fv3jedi_field_registry%get(c_key_out,xout)

call fv3jedi_varchange_multiplyadjoint(self,xin,xout)

end subroutine c_fv3jedi_varchange_multiplyadjoint

! ----------------------------------------------------------------------------

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_varchange_multiplyinverse(c_key_self, c_key_in, c_key_out) &
           bind (c,name='fv3jedi_varchange_multiplyinverse_f90')

use iso_c_binding
use fv3jedi_varchange_mod
use fv3jedi_fields_mod, only: fv3jedi_field_registry, fv3jedi_field
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_varchange), pointer :: self
type(fv3jedi_field), pointer :: xin
type(fv3jedi_field), pointer :: xout

call fv3jedi_varchange_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_in,xin)
call fv3jedi_field_registry%get(c_key_out,xout)

call fv3jedi_varchange_multiplyinverse(self,xin,xout)

end subroutine c_fv3jedi_varchange_multiplyinverse

! ----------------------------------------------------------------------------

subroutine c_fv3jedi_varchange_multiplyinverseadjoint(c_key_self, c_key_in, &
      c_key_out) bind (c,name='fv3jedi_varchange_multiplyinverseadjoint_f90')

use iso_c_binding
use fv3jedi_varchange_mod
use fv3jedi_fields_mod, only: fv3jedi_field_registry, fv3jedi_field
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_varchange), pointer :: self
type(fv3jedi_field), pointer :: xin
type(fv3jedi_field), pointer :: xout

call fv3jedi_varchange_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_in,xin)
call fv3jedi_field_registry%get(c_key_out,xout)

call fv3jedi_varchange_multiplyinverseadjoint(self,xin,xout)

end subroutine c_fv3jedi_varchange_multiplyinverseadjoint

! ----------------------------------------------------------------------------
