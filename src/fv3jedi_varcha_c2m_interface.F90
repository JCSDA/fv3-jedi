! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_varcha_c2m_setup(c_key_self, c_key_state_bg, c_key_state_fg, &
           c_key_geom, c_conf) bind (c,name='fv3jedi_varcha_c2m_setup_f90')

use iso_c_binding
use fv3jedi_varcha_c2m_mod
use fv3jedi_fields_mod
use fv3jedi_geom_mod

implicit none
integer(c_int), intent(inout) :: c_key_self     !< Change variable structure
integer(c_int), intent(in)    :: c_key_state_bg !< Background key
integer(c_int), intent(in)    :: c_key_state_fg !< First guess key
integer(c_int), intent(in)    :: c_key_geom     !< Geom key
type(c_ptr),    intent(in)    :: c_conf         !< Configuration

type(fv3jedi_varcha_c2m), pointer :: self
type(fv3jedi_field), pointer :: bg
type(fv3jedi_field), pointer :: fg
type(fv3jedi_geom), pointer :: geom

call fv3jedi_varcha_c2m_registry%init()
call fv3jedi_varcha_c2m_registry%add(c_key_self)
call fv3jedi_varcha_c2m_registry%get(c_key_self, self)

call fv3jedi_field_registry%get(c_key_state_bg,bg)
call fv3jedi_field_registry%get(c_key_state_fg,fg)

call fv3jedi_geom_registry%get(c_key_geom,geom)

call fv3jedi_varcha_c2m_setup(self, bg, fg, geom, c_conf)

end subroutine c_fv3jedi_varcha_c2m_setup

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_varcha_c2m_delete(c_key_self) &
           bind (c,name='fv3jedi_varcha_c2m_delete_f90')

use iso_c_binding
use fv3jedi_varcha_c2m_mod

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Change variable structure

type(fv3jedi_varcha_c2m), pointer :: self

call fv3jedi_varcha_c2m_registry%get(c_key_self,self)
call fv3jedi_varcha_c2m_delete(self)
call fv3jedi_varcha_c2m_registry%remove(c_key_self)

end subroutine c_fv3jedi_varcha_c2m_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_varcha_c2m_multiply(c_key_self, c_key_in, c_key_out) &
           bind (c,name='fv3jedi_varcha_c2m_multiply_f90')

use iso_c_binding
use fv3jedi_varcha_c2m_mod
use fv3jedi_fields_mod, only: fv3jedi_field_registry, fv3jedi_field
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_varcha_c2m), pointer :: self
type(fv3jedi_field), pointer :: xin
type(fv3jedi_field), pointer :: xout

call fv3jedi_varcha_c2m_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_in,xin)
call fv3jedi_field_registry%get(c_key_out,xout)

call fv3jedi_varcha_c2m_multiply(self,xin,xout)

end subroutine c_fv3jedi_varcha_c2m_multiply

! ----------------------------------------------------------------------------

subroutine c_fv3jedi_varcha_c2m_multiplyadjoint(c_key_self, c_key_in, &
           c_key_out) bind (c,name='fv3jedi_varcha_c2m_multiplyadjoint_f90')

use iso_c_binding
use fv3jedi_varcha_c2m_mod
use fv3jedi_fields_mod, only: fv3jedi_field_registry, fv3jedi_field
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_varcha_c2m), pointer :: self
type(fv3jedi_field), pointer :: xin
type(fv3jedi_field), pointer :: xout

call fv3jedi_varcha_c2m_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_in,xin)
call fv3jedi_field_registry%get(c_key_out,xout)

call fv3jedi_varcha_c2m_multiplyadjoint(self,xin,xout)

end subroutine c_fv3jedi_varcha_c2m_multiplyadjoint

! ----------------------------------------------------------------------------

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_varcha_c2m_multiplyinverse(c_key_self, c_key_in, c_key_out) &
           bind (c,name='fv3jedi_varcha_c2m_multiplyinverse_f90')

use iso_c_binding
use fv3jedi_varcha_c2m_mod
use fv3jedi_fields_mod, only: fv3jedi_field_registry, fv3jedi_field
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_varcha_c2m), pointer :: self
type(fv3jedi_field), pointer :: xin
type(fv3jedi_field), pointer :: xout

call fv3jedi_varcha_c2m_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_in,xin)
call fv3jedi_field_registry%get(c_key_out,xout)

call fv3jedi_varcha_c2m_multiplyinverse(self,xin,xout)

end subroutine c_fv3jedi_varcha_c2m_multiplyinverse

! ----------------------------------------------------------------------------

subroutine c_fv3jedi_varcha_c2m_multiplyinverseadjoint(c_key_self, c_key_in, &
      c_key_out) bind (c,name='fv3jedi_varcha_c2m_multiplyinverseadjoint_f90')

use iso_c_binding
use fv3jedi_varcha_c2m_mod
use fv3jedi_fields_mod, only: fv3jedi_field_registry, fv3jedi_field
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out

type(fv3jedi_varcha_c2m), pointer :: self
type(fv3jedi_field), pointer :: xin
type(fv3jedi_field), pointer :: xout

call fv3jedi_varcha_c2m_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_in,xin)
call fv3jedi_field_registry%get(c_key_out,xout)

call fv3jedi_varcha_c2m_multiplyinverseadjoint(self,xin,xout)

end subroutine c_fv3jedi_varcha_c2m_multiplyinverseadjoint

! ----------------------------------------------------------------------------
