! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_create_c(c_key_self, c_key_geom, c_vars) bind(c,name='fv3jedi_field_create_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use fv3jedi_geom_mod
use ufo_vars_mod
implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in) :: c_key_geom !< Geometry
type(c_ptr), intent(in)    :: c_vars     !< List of variables

type(fv3jedi_field), pointer :: self
type(fv3jedi_geom),  pointer :: geom
type(ufo_vars) :: vars

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_field_registry%init()
call fv3jedi_field_registry%add(c_key_self)
call fv3jedi_field_registry%get(c_key_self,self)

call create(self, geom, vars)

end subroutine fv3jedi_field_create_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_delete_c(c_key_self) bind(c,name='fv3jedi_field_delete_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
implicit none
integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_field), pointer :: self

call fv3jedi_field_registry%get(c_key_self,self)

call delete(self)

call fv3jedi_field_registry%remove(c_key_self)

end subroutine fv3jedi_field_delete_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_zero_c(c_key_self) bind(c,name='fv3jedi_field_zero_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
implicit none
integer(c_int), intent(in) :: c_key_self
type(fv3jedi_field), pointer :: self

call fv3jedi_field_registry%get(c_key_self,self)
call zeros(self)

end subroutine fv3jedi_field_zero_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_dirac_c(c_key_self,c_conf) bind(c,name='fv3jedi_field_dirac_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in)    :: c_conf !< Configuration
type(fv3jedi_field), pointer :: self

call fv3jedi_field_registry%get(c_key_self,self)
call dirac(self,c_conf)

end subroutine fv3jedi_field_dirac_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_random_c(c_key_self) bind(c,name='fv3jedi_field_random_f90')
use iso_c_binding
use fv3jedi_fields_mod, only: fv3jedi_field_registry, random
use fv3jedi_fields_utils_mod, only: fv3jedi_field
implicit none
integer(c_int), intent(in) :: c_key_self
type(fv3jedi_field), pointer :: self

call fv3jedi_field_registry%get(c_key_self,self)
call random(self)

end subroutine fv3jedi_field_random_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_convert_to_c(c_key_fld, c_key_ug) bind (c,name='fv3jedi_field_convert_to_f90')

use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use unstructured_grid_mod
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(c_int), intent(in) :: c_key_ug
type(fv3jedi_field), pointer :: fld
type(unstructured_grid), pointer :: ug

call fv3jedi_field_registry%get(c_key_fld,fld)
call unstructured_grid_registry%get(c_key_ug,ug)

call convert_to_ug(fld, ug)

end subroutine fv3jedi_field_convert_to_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_convert_from_c(c_key_fld, c_key_ug) bind (c,name='fv3jedi_field_convert_from_f90')

use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use unstructured_grid_mod
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(c_int), intent(in) :: c_key_ug
type(fv3jedi_field), pointer :: fld
type(unstructured_grid), pointer :: ug

call fv3jedi_field_registry%get(c_key_fld,fld)
call unstructured_grid_registry%get(c_key_ug,ug)

call convert_from_ug(fld, ug)

end subroutine fv3jedi_field_convert_from_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_copy_c(c_key_self,c_key_rhs) bind(c,name='fv3jedi_field_copy_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_field), pointer :: self
type(fv3jedi_field), pointer :: rhs
call fv3jedi_field_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_rhs,rhs)

call copy(self, rhs)

end subroutine fv3jedi_field_copy_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_self_add_c(c_key_self,c_key_rhs) bind(c,name='fv3jedi_field_self_add_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_field), pointer :: self
type(fv3jedi_field), pointer :: rhs
call fv3jedi_field_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_rhs,rhs)

call self_add(self,rhs)

end subroutine fv3jedi_field_self_add_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_self_schur_c(c_key_self,c_key_rhs) bind(c,name='fv3jedi_field_self_schur_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_field), pointer :: self
type(fv3jedi_field), pointer :: rhs
call fv3jedi_field_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_rhs,rhs)

call self_schur(self,rhs)

end subroutine fv3jedi_field_self_schur_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_self_sub_c(c_key_self,c_key_rhs) bind(c,name='fv3jedi_field_self_sub_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_field), pointer :: self
type(fv3jedi_field), pointer :: rhs
call fv3jedi_field_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_rhs,rhs)

call self_sub(self,rhs)

end subroutine fv3jedi_field_self_sub_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_self_mul_c(c_key_self,c_zz) bind(c,name='fv3jedi_field_self_mul_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use kinds
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
type(fv3jedi_field), pointer :: self
real(kind=kind_real) :: zz

call fv3jedi_field_registry%get(c_key_self,self)
zz = c_zz

call self_mul(self,zz)

end subroutine fv3jedi_field_self_mul_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_axpy_c(c_key_self,c_zz,c_key_rhs) bind(c,name='fv3jedi_field_axpy_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use kinds
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_field), pointer :: self
type(fv3jedi_field), pointer :: rhs
real(kind=kind_real) :: zz

call fv3jedi_field_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_rhs,rhs)
zz = c_zz

call axpy(self,zz,rhs)

end subroutine fv3jedi_field_axpy_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_dot_prod_c(c_key_fld1,c_key_fld2,c_prod) bind(c,name='fv3jedi_field_dot_prod_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use kinds
implicit none
integer(c_int), intent(in)    :: c_key_fld1, c_key_fld2
real(c_double), intent(inout) :: c_prod
real(kind=kind_real) :: zz
type(fv3jedi_field), pointer :: fld1, fld2

call fv3jedi_field_registry%get(c_key_fld1,fld1)
call fv3jedi_field_registry%get(c_key_fld2,fld2)

call dot_prod(fld1,fld2,zz)

c_prod = zz

end subroutine fv3jedi_field_dot_prod_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_add_incr_c(c_key_self,c_key_rhs) bind(c,name='fv3jedi_field_add_incr_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs
type(fv3jedi_field), pointer :: self
type(fv3jedi_field), pointer :: rhs

call fv3jedi_field_registry%get(c_key_self,self)
call fv3jedi_field_registry%get(c_key_rhs,rhs)

call add_incr(self,rhs)

end subroutine fv3jedi_field_add_incr_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_diff_incr_c(c_key_lhs,c_key_x1,c_key_x2) bind(c,name='fv3jedi_field_diff_incr_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
implicit none
integer(c_int), intent(in) :: c_key_lhs
integer(c_int), intent(in) :: c_key_x1
integer(c_int), intent(in) :: c_key_x2
type(fv3jedi_field), pointer :: lhs
type(fv3jedi_field), pointer :: x1
type(fv3jedi_field), pointer :: x2

call fv3jedi_field_registry%get(c_key_lhs,lhs)
call fv3jedi_field_registry%get(c_key_x1,x1)
call fv3jedi_field_registry%get(c_key_x2,x2)

call diff_incr(lhs,x1,x2)

end subroutine fv3jedi_field_diff_incr_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_change_resol_c(c_key_fld,c_key_rhs) bind(c,name='fv3jedi_field_change_resol_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(c_int), intent(in) :: c_key_rhs
type(fv3jedi_field), pointer :: fld, rhs

call fv3jedi_field_registry%get(c_key_fld,fld)
call fv3jedi_field_registry%get(c_key_rhs,rhs)

call change_resol(fld,rhs)

end subroutine fv3jedi_field_change_resol_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_read_file_c(c_key_fld, c_conf, c_dt) bind(c,name='fv3jedi_field_read_file_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use datetime_mod

implicit none
integer(c_int), intent(in) :: c_key_fld  !< Fields
type(c_ptr), intent(in)    :: c_conf !< Configuration
type(c_ptr), intent(inout) :: c_dt   !< DateTime

type(fv3jedi_field), pointer :: fld
type(datetime) :: fdate

call fv3jedi_field_registry%get(c_key_fld,fld)
call c_f_datetime(c_dt, fdate)
call read_file(fld, c_conf, fdate)

end subroutine fv3jedi_field_read_file_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_analytic_init_c(c_key_fld, c_key_geom, c_conf, c_dt) bind(c,name='fv3jedi_field_analytic_init_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use fv3jedi_geom_mod
use datetime_mod

implicit none
integer(c_int), intent(in) :: c_key_fld  !< Fields
integer(c_int), intent(in) :: c_key_geom  !< Geometry
type(c_ptr), intent(in)    :: c_conf !< Configuration
type(c_ptr), intent(inout) :: c_dt   !< DateTime

type(fv3jedi_field), pointer :: fld
type(fv3jedi_geom), pointer :: geom
type(datetime) :: fdate

call fv3jedi_field_registry%get(c_key_fld,fld)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call c_f_datetime(c_dt, fdate)
call analytic_IC(fld, geom, c_conf, fdate)

end subroutine fv3jedi_field_analytic_init_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_write_file_c(c_key_fld, c_conf, c_dt) bind(c,name='fv3jedi_field_write_file_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use datetime_mod

implicit none
integer(c_int), intent(in) :: c_key_fld  !< Fields
type(c_ptr), intent(in) :: c_conf !< Configuration
type(c_ptr), intent(in) :: c_dt   !< DateTime

type(fv3jedi_field), pointer :: fld
type(datetime) :: fdate

call fv3jedi_field_registry%get(c_key_fld,fld)
call c_f_datetime(c_dt, fdate)
call write_file(fld, c_conf, fdate)

end subroutine fv3jedi_field_write_file_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_gpnorm_c(c_key_fld, kf, pstat) bind(c,name='fv3jedi_field_gpnorm_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use kinds
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(c_int), intent(in) :: kf
real(c_double), intent(inout) :: pstat(3*kf)

type(fv3jedi_field), pointer :: fld
real(kind=kind_real) :: zstat(3, kf)
integer :: jj, js, jf

call fv3jedi_field_registry%get(c_key_fld,fld)

call gpnorm(fld, kf, zstat)
jj=0
do jf = 1, kf
  do js = 1, 3
    jj=jj+1
    pstat(jj) = zstat(js,jf)
  enddo
enddo

end subroutine fv3jedi_field_gpnorm_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_rms_c(c_key_fld, prms) bind(c,name='fv3jedi_field_rms_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use kinds
implicit none
integer(c_int), intent(in) :: c_key_fld
real(c_double), intent(inout) :: prms

type(fv3jedi_field), pointer :: fld
real(kind=kind_real) :: zz

call fv3jedi_field_registry%get(c_key_fld,fld)

call fldrms(fld, zz)

prms = zz

end subroutine fv3jedi_field_rms_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_interp_c(c_key_fld,c_key_loc,c_vars,c_key_gom) bind(c,name='fv3jedi_field_interp_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use ioda_locs_mod
use ioda_locs_mod_c, only: ioda_locs_registry
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
implicit none
integer(c_int), intent(in) :: c_key_fld  !< Fields to be interpolated
integer(c_int), intent(in) :: c_key_loc  !< List of requested locations
type(c_ptr), intent(in)    :: c_vars     !< List of requested variables
integer(c_int), intent(in) :: c_key_gom  !< Interpolated values
type(fv3jedi_field), pointer :: fld
type(ioda_locs),  pointer :: locs
type(ufo_geovals),  pointer :: gom
type(ufo_vars) :: vars

call ufo_vars_setup(vars, c_vars)

call fv3jedi_field_registry%get(c_key_fld, fld)
call ioda_locs_registry%get(c_key_loc, locs)
call ufo_geovals_registry%get(c_key_gom, gom)

call interp(fld, locs, vars, gom)

end subroutine fv3jedi_field_interp_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_interp_tl_c(c_key_fld,c_key_loc,c_vars,c_key_gom) bind(c,name='fv3jedi_field_interp_tl_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use ioda_locs_mod
use ioda_locs_mod_c, only: ioda_locs_registry
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
implicit none
integer(c_int), intent(in) :: c_key_fld  !< Fields to be interpolated
integer(c_int), intent(in) :: c_key_loc  !< List of requested locations
type(c_ptr), intent(in)    :: c_vars     !< List of requested variables
integer(c_int), intent(in) :: c_key_gom  !< Interpolated values
type(fv3jedi_field), pointer :: fld
type(ioda_locs),  pointer :: locs
type(ufo_geovals),  pointer :: gom
type(ufo_vars) :: vars

call ufo_vars_setup(vars, c_vars)

call fv3jedi_field_registry%get(c_key_fld, fld)
call ioda_locs_registry%get(c_key_loc, locs)
call ufo_geovals_registry%get(c_key_gom, gom)

call interp_tl(fld, locs, vars, gom)

end subroutine fv3jedi_field_interp_tl_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_interp_ad_c(c_key_fld,c_key_loc,c_vars,c_key_gom) bind(c,name='fv3jedi_field_interp_ad_f90')
use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field
use ioda_locs_mod
use ioda_locs_mod_c, only: ioda_locs_registry
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
implicit none
integer(c_int), intent(in) :: c_key_fld  !< Fields to be interpolated
integer(c_int), intent(in) :: c_key_loc  !< List of requested locations
type(c_ptr), intent(in)    :: c_vars     !< List of requested variables
integer(c_int), intent(in) :: c_key_gom  !< Interpolated values
type(fv3jedi_field), pointer :: fld
type(ioda_locs),  pointer :: locs
type(ufo_geovals),  pointer :: gom
type(ufo_vars) :: vars

call ufo_vars_setup(vars, c_vars)

call fv3jedi_field_registry%get(c_key_fld, fld)
call ioda_locs_registry%get(c_key_loc, locs)
call ufo_geovals_registry%get(c_key_gom, gom)

call interp_ad(fld, locs, vars, gom)

end subroutine fv3jedi_field_interp_ad_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_field_sizes_c(c_key_self,nx,ny,nf) bind(c,name='fv3jedi_field_sizes_f90')

use iso_c_binding
use fv3jedi_fields_mod
use fv3jedi_fields_utils_mod, only: fv3jedi_field

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: nx,ny,nf
type(fv3jedi_field), pointer :: self

call fv3jedi_field_registry%get(c_key_self,self)

nf = self%nf
nx = self%geom%npx
ny = self%geom%npy

end subroutine fv3jedi_field_sizes_c

! ------------------------------------------------------------------------------   
