! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_b_setup(c_key_self, c_conf, c_key_geom) &
          & bind (c,name='fv3jedi_b_setup_f90')

use iso_c_binding
use fv3jedi_covariance_mod
use fv3jedi_geom_mod, only: fv3jedi_geom, fv3jedi_geom_registry

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Background error covariance structure
type(c_ptr), intent(in)    :: c_conf         !< Configuration
integer(c_int), intent(in) :: c_key_geom     !< Geometry
type(fv3jedi_covar), pointer :: self
type(fv3jedi_geom),  pointer :: geom

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_covar_registry%init()
call fv3jedi_covar_registry%add(c_key_self)
call fv3jedi_covar_registry%get(c_key_self, self)

call fv3jedi_covar_setup(self, geom, c_conf)

end subroutine c_fv3jedi_b_setup

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_b_delete(c_key_self) bind (c,name='fv3jedi_b_delete_f90')

use iso_c_binding
use fv3jedi_covariance_mod

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Background error covariance structure
type(fv3jedi_covar), pointer :: self

call fv3jedi_covar_registry%get(c_key_self,self)
call fv3jedi_covar_delete(self)
call fv3jedi_covar_registry%remove(c_key_self)

end subroutine c_fv3jedi_b_delete

! ------------------------------------------------------------------------------

!> Multiply streamfunction by inverse of covariance

subroutine c_fv3jedi_b_inv_mult(c_key_self, c_key_in, c_key_out) bind(c,name='fv3jedi_b_invmult_f90')

use iso_c_binding
use fv3jedi_covariance_mod
use fv3jedi_fields_mod, only: fv3jedi_field_registry
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out
type(fv3jedi_covar), pointer :: self
type(fv3jedi_field), pointer :: xin
type(fv3jedi_field), pointer :: xout

call fv3jedi_covar_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_in,xin)
call fv3jedi_field_registry%get(c_key_out,xout)

!call fv3jedi_covar_sqrt_inv_mult(self%nx,self%ny,xctl,xin,self)
!call zeros(xout)
!call fv3jedi_covar_sqrt_inv_mult_ad(self%nx,self%ny,xctl,xout,self)

end subroutine c_fv3jedi_b_inv_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by covariance

subroutine c_fv3jedi_b_mult(c_key_self, c_key_in, c_key_out) bind(c,name='fv3jedi_b_mult_f90')

use iso_c_binding
use fv3jedi_covariance_mod
use fv3jedi_fields_mod, only: fv3jedi_field_registry
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_in
integer(c_int), intent(in) :: c_key_out
type(fv3jedi_covar), pointer :: self
type(fv3jedi_field), pointer :: xin
type(fv3jedi_field), pointer :: xout

call fv3jedi_covar_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_in,xin)
call fv3jedi_field_registry%get(c_key_out,xout)

!call fv3jedi_covar_sqrt_mult_ad(self%nx,self%ny,xin,xctl,self)
!call zeros(xout)
!call fv3jedi_covar_sqrt_mult(self%nx,self%ny,xout,xctl,self)

end subroutine c_fv3jedi_b_mult

! ------------------------------------------------------------------------------

!> Generate randomized increment

subroutine c_fv3jedi_b_randomize(c_key_self, c_key_out) bind(c,name='fv3jedi_b_randomize_f90')

use iso_c_binding
use fv3jedi_covariance_mod
use fv3jedi_fields_mod, only: fv3jedi_field_registry, random
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use kinds

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_out
type(fv3jedi_covar), pointer :: self
type(fv3jedi_field), pointer :: xout

call fv3jedi_covar_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_out,xout)

call random(xout)

end subroutine c_fv3jedi_b_randomize

! ------------------------------------------------------------------------------
