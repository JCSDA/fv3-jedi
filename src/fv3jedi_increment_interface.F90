! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_create_c(c_key_self, c_key_geom, c_vars) bind(c,name='fv3jedi_increment_create_f90')
use iso_c_binding
use fv3jedi_increment_mod
use fv3jedi_geom_mod
use fv3jedi_vars_mod

implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in) :: c_key_geom !< Geometry
type(c_ptr), intent(in)    :: c_vars     !< List of variables

type(fv3jedi_increment), pointer :: self
type(fv3jedi_geom),  pointer :: geom
type(fv3jedi_vars) :: vars

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_increment_registry%init()
call fv3jedi_increment_registry%add(c_key_self)
call fv3jedi_increment_registry%get(c_key_self,self)

call fv3jedi_vars_create(c_vars,vars)
call create(self, geom, vars)

end subroutine fv3jedi_increment_create_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_delete_c(c_key_self) bind(c,name='fv3jedi_increment_delete_f90')
use iso_c_binding
use fv3jedi_increment_mod
implicit none
integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_increment), pointer :: self

call fv3jedi_increment_registry%get(c_key_self,self)

call delete(self)

call fv3jedi_increment_registry%remove(c_key_self)

end subroutine fv3jedi_increment_delete_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_zero_c(c_key_self) bind(c,name='fv3jedi_increment_zero_f90')
use iso_c_binding
use fv3jedi_increment_mod
implicit none
integer(c_int), intent(in) :: c_key_self
type(fv3jedi_increment), pointer :: self

call fv3jedi_increment_registry%get(c_key_self,self)
call zeros(self)

end subroutine fv3jedi_increment_zero_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_dirac_c(c_key_self,c_conf,c_key_geom) bind(c,name='fv3jedi_increment_dirac_f90')
use iso_c_binding
use fv3jedi_increment_mod
use fv3jedi_geom_mod
implicit none
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in)    :: c_conf !< Configuration
integer(c_int), intent(in) :: c_key_geom !< Geometry
type(fv3jedi_increment), pointer :: self
type(fv3jedi_geom),  pointer :: geom

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_increment_registry%get(c_key_self,self)
call dirac(self,c_conf,geom)

end subroutine fv3jedi_increment_dirac_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_random_c(c_key_self) bind(c,name='fv3jedi_increment_random_f90')
use iso_c_binding
use fv3jedi_increment_mod
implicit none
integer(c_int), intent(in) :: c_key_self
type(fv3jedi_increment), pointer :: self

call fv3jedi_increment_registry%get(c_key_self,self)
call random(self)

end subroutine fv3jedi_increment_random_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_ug_coord_c(c_key_inc, c_key_ug, c_colocated, c_key_geom) bind (c,name='fv3jedi_increment_ug_coord_f90')

use iso_c_binding
use fv3jedi_increment_mod
use unstructured_grid_mod
use fv3jedi_geom_mod
implicit none
integer(c_int), intent(in) :: c_key_inc
integer(c_int), intent(in) :: c_key_ug
integer(c_int), intent(in) :: c_colocated
integer(c_int), intent(in) :: c_key_geom !< Geometry
type(fv3jedi_increment), pointer :: inc
type(unstructured_grid), pointer :: ug
integer :: colocated
type(fv3jedi_geom),  pointer :: geom

colocated = c_colocated

call fv3jedi_increment_registry%get(c_key_inc,inc)
call unstructured_grid_registry%get(c_key_ug,ug)
call fv3jedi_geom_registry%get(c_key_geom, geom)

call ug_coord(inc, ug, colocated, geom)

end subroutine fv3jedi_increment_ug_coord_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_increment_to_ug_c(c_key_inc, c_key_ug, c_colocated) bind (c,name='fv3jedi_increment_increment_to_ug_f90')

use iso_c_binding
use fv3jedi_increment_mod
use unstructured_grid_mod
implicit none
integer(c_int), intent(in) :: c_key_inc
integer(c_int), intent(in) :: c_key_ug
integer(c_int), intent(in) :: c_colocated
type(fv3jedi_increment), pointer :: inc
type(unstructured_grid), pointer :: ug
integer :: colocated

colocated = c_colocated

call fv3jedi_increment_registry%get(c_key_inc,inc)
call unstructured_grid_registry%get(c_key_ug,ug)

call increment_to_ug(inc, ug, colocated)

end subroutine fv3jedi_increment_increment_to_ug_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_increment_from_ug_c(c_key_inc, c_key_ug) bind (c,name='fv3jedi_increment_increment_from_ug_f90')

use iso_c_binding
use fv3jedi_increment_mod
use unstructured_grid_mod
implicit none
integer(c_int), intent(in) :: c_key_inc
integer(c_int), intent(in) :: c_key_ug
type(fv3jedi_increment), pointer :: inc
type(unstructured_grid), pointer :: ug

call fv3jedi_increment_registry%get(c_key_inc,inc)
call unstructured_grid_registry%get(c_key_ug,ug)

call increment_from_ug(inc, ug)

end subroutine fv3jedi_increment_increment_from_ug_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_copy_c(c_key_self,c_key_rhs) bind(c,name='fv3jedi_increment_copy_f90')
use iso_c_binding
use fv3jedi_increment_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_increment), pointer :: self
type(fv3jedi_increment), pointer :: rhs
call fv3jedi_increment_registry%get(c_key_self,self)
call fv3jedi_increment_registry%get(c_key_rhs,rhs)

call copy(self, rhs)

end subroutine fv3jedi_increment_copy_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_self_add_c(c_key_self,c_key_rhs) bind(c,name='fv3jedi_increment_self_add_f90')
use iso_c_binding
use fv3jedi_increment_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_increment), pointer :: self
type(fv3jedi_increment), pointer :: rhs
call fv3jedi_increment_registry%get(c_key_self,self)
call fv3jedi_increment_registry%get(c_key_rhs,rhs)

call self_add(self,rhs)

end subroutine fv3jedi_increment_self_add_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_self_schur_c(c_key_self,c_key_rhs) bind(c,name='fv3jedi_increment_self_schur_f90')
use iso_c_binding
use fv3jedi_increment_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_increment), pointer :: self
type(fv3jedi_increment), pointer :: rhs
call fv3jedi_increment_registry%get(c_key_self,self)
call fv3jedi_increment_registry%get(c_key_rhs,rhs)

call self_schur(self,rhs)

end subroutine fv3jedi_increment_self_schur_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_self_sub_c(c_key_self,c_key_rhs) bind(c,name='fv3jedi_increment_self_sub_f90')
use iso_c_binding
use fv3jedi_increment_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_increment), pointer :: self
type(fv3jedi_increment), pointer :: rhs
call fv3jedi_increment_registry%get(c_key_self,self)
call fv3jedi_increment_registry%get(c_key_rhs,rhs)

call self_sub(self,rhs)

end subroutine fv3jedi_increment_self_sub_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_self_mul_c(c_key_self,c_zz) bind(c,name='fv3jedi_increment_self_mul_f90')
use iso_c_binding
use fv3jedi_increment_mod
use kinds
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
type(fv3jedi_increment), pointer :: self
real(kind=kind_real) :: zz

call fv3jedi_increment_registry%get(c_key_self,self)
zz = c_zz

call self_mul(self,zz)

end subroutine fv3jedi_increment_self_mul_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_axpy_inc_c(c_key_self,c_zz,c_key_rhs) bind(c,name='fv3jedi_increment_axpy_inc_f90')
use iso_c_binding
use fv3jedi_increment_mod
use kinds
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_increment), pointer :: self
type(fv3jedi_increment), pointer :: rhs
real(kind=kind_real) :: zz

call fv3jedi_increment_registry%get(c_key_self,self)
call fv3jedi_increment_registry%get(c_key_rhs,rhs)
zz = c_zz

call axpy_inc(self,zz,rhs)

end subroutine fv3jedi_increment_axpy_inc_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_axpy_state_c(c_key_self,c_zz,c_key_rhs) bind(c,name='fv3jedi_increment_axpy_state_f90')
use iso_c_binding
use fv3jedi_increment_mod
use fv3jedi_state_mod, only: fv3jedi_state, fv3jedi_state_registry
use kinds
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_increment), pointer :: self
type(fv3jedi_state), pointer :: rhs
real(kind=kind_real) :: zz

call fv3jedi_increment_registry%get(c_key_self,self)
call fv3jedi_state_registry%get(c_key_rhs,rhs)
zz = c_zz

call axpy_state(self,zz,rhs)

end subroutine fv3jedi_increment_axpy_state_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_dot_prod_c(c_key_inc1,c_key_inc2,c_prod) bind(c,name='fv3jedi_increment_dot_prod_f90')
use iso_c_binding
use fv3jedi_increment_mod
use kinds
implicit none
integer(c_int), intent(in)    :: c_key_inc1, c_key_inc2
real(c_double), intent(inout) :: c_prod
real(kind=kind_real) :: zz
type(fv3jedi_increment), pointer :: inc1, inc2

call fv3jedi_increment_registry%get(c_key_inc1,inc1)
call fv3jedi_increment_registry%get(c_key_inc2,inc2)

call dot_prod(inc1,inc2,zz)

c_prod = zz

end subroutine fv3jedi_increment_dot_prod_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_add_incr_c(c_key_self,c_key_rhs) bind(c,name='fv3jedi_increment_add_incr_f90')
use iso_c_binding
use fv3jedi_increment_mod
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs
type(fv3jedi_increment), pointer :: self
type(fv3jedi_increment), pointer :: rhs

call fv3jedi_increment_registry%get(c_key_self,self)
call fv3jedi_increment_registry%get(c_key_rhs,rhs)

call add_incr(self,rhs)

end subroutine fv3jedi_increment_add_incr_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_diff_incr_c(c_key_lhs,c_key_x1,c_key_x2) bind(c,name='fv3jedi_increment_diff_incr_f90')
use iso_c_binding
use fv3jedi_increment_mod
use fv3jedi_state_mod, only: fv3jedi_state, fv3jedi_state_registry
implicit none
integer(c_int), intent(in) :: c_key_lhs
integer(c_int), intent(in) :: c_key_x1
integer(c_int), intent(in) :: c_key_x2
type(fv3jedi_increment), pointer :: lhs
type(fv3jedi_state), pointer :: x1
type(fv3jedi_state), pointer :: x2

call fv3jedi_increment_registry%get(c_key_lhs,lhs)
call fv3jedi_state_registry%get(c_key_x1,x1)
call fv3jedi_state_registry%get(c_key_x2,x2)

call diff_incr(lhs,x1,x2)

end subroutine fv3jedi_increment_diff_incr_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_change_resol_c(c_key_inc,c_key_rhs) bind(c,name='fv3jedi_increment_change_resol_f90')
use iso_c_binding
use fv3jedi_increment_mod
implicit none
integer(c_int), intent(in) :: c_key_inc
integer(c_int), intent(in) :: c_key_rhs
type(fv3jedi_increment), pointer :: inc, rhs

call fv3jedi_increment_registry%get(c_key_inc,inc)
call fv3jedi_increment_registry%get(c_key_rhs,rhs)

call change_resol(inc,rhs)

end subroutine fv3jedi_increment_change_resol_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_read_file_c(c_key_geom, c_key_inc, c_conf, c_dt) bind(c,name='fv3jedi_increment_read_file_f90')
use iso_c_binding
use fv3jedi_increment_mod
use datetime_mod
use fv3jedi_geom_mod

implicit none
integer(c_int), intent(in) :: c_key_inc  !< Increment
type(c_ptr), intent(in)    :: c_conf !< Configuration
type(c_ptr), intent(inout) :: c_dt   !< DateTime
integer(c_int), intent(in) :: c_key_geom  !< Geometry

type(fv3jedi_increment), pointer :: inc
type(datetime) :: fdate
type(fv3jedi_geom),  pointer :: geom

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_increment_registry%get(c_key_inc,inc)
call c_f_datetime(c_dt, fdate)
call read_file(geom, inc, c_conf, fdate)

end subroutine fv3jedi_increment_read_file_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_write_file_c(c_key_geom, c_key_inc, c_conf, c_dt) bind(c,name='fv3jedi_increment_write_file_f90')
use iso_c_binding
use fv3jedi_increment_mod
use datetime_mod
use fv3jedi_geom_mod

implicit none
integer(c_int), intent(in) :: c_key_inc  !< Increment
type(c_ptr), intent(in) :: c_conf !< Configuration
type(c_ptr), intent(in) :: c_dt   !< DateTime
integer(c_int), intent(in) :: c_key_geom  !< Geometry

type(fv3jedi_increment), pointer :: inc
type(datetime) :: fdate
type(fv3jedi_geom),  pointer :: geom

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_increment_registry%get(c_key_inc,inc)
call c_f_datetime(c_dt, fdate)
call write_file(geom, inc, c_conf, fdate)

end subroutine fv3jedi_increment_write_file_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_gpnorm_c(c_key_inc, kf, pstat) bind(c,name='fv3jedi_increment_gpnorm_f90')
use iso_c_binding
use fv3jedi_increment_mod
use kinds
implicit none
integer(c_int), intent(in) :: c_key_inc
integer(c_int), intent(in) :: kf
real(c_double), intent(inout) :: pstat(3*kf)

type(fv3jedi_increment), pointer :: inc
real(kind=kind_real) :: zstat(3, kf)
integer :: jj, js, jf

call fv3jedi_increment_registry%get(c_key_inc,inc)

call gpnorm(inc, kf, zstat)
jj=0
do jf = 1, kf
  do js = 1, 3
    jj=jj+1
    pstat(jj) = zstat(js,jf)
  enddo
enddo

end subroutine fv3jedi_increment_gpnorm_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_rms_c(c_key_inc, prms) bind(c,name='fv3jedi_increment_rms_f90')
use iso_c_binding
use fv3jedi_increment_mod
use kinds
implicit none
integer(c_int), intent(in) :: c_key_inc
real(c_double), intent(inout) :: prms

type(fv3jedi_increment), pointer :: inc
real(kind=kind_real) :: zz

call fv3jedi_increment_registry%get(c_key_inc,inc)

call incrms(inc, zz)

prms = zz

end subroutine fv3jedi_increment_rms_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_getvalues_tl_c(c_key_geom, c_key_inc,c_key_loc,c_vars,c_key_gom,c_key_traj) bind(c,name='fv3jedi_increment_getvalues_tl_f90')
use iso_c_binding
use fv3jedi_increment_mod
use ioda_locs_mod
use ioda_locs_mod_c, only: ioda_locs_registry
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
use fv3jedi_getvaltraj_mod, only: fv3jedi_getvaltraj, fv3jedi_getvaltraj_registry
use fv3jedi_geom_mod
implicit none
integer(c_int), intent(in) :: c_key_inc  !< Increment to be interpolated
integer(c_int), intent(in) :: c_key_loc  !< List of requested locations
type(c_ptr), intent(in)    :: c_vars     !< List of requested variables
integer(c_int), intent(in) :: c_key_gom  !< Interpolated values
integer(c_int), intent(in) :: c_key_traj !< Trajectory for interpolation/transforms
integer(c_int), intent(in) :: c_key_geom  !< Geometry
type(fv3jedi_increment), pointer :: inc
type(ioda_locs),  pointer :: locs
type(ufo_geovals),  pointer :: gom
type(ufo_vars) :: vars
type(fv3jedi_getvaltraj), pointer :: traj
type(fv3jedi_geom),  pointer :: geom

call ufo_vars_setup(vars, c_vars)

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_increment_registry%get(c_key_inc, inc)
call ioda_locs_registry%get(c_key_loc, locs)
call ufo_geovals_registry%get(c_key_gom, gom)
call fv3jedi_getvaltraj_registry%get(c_key_traj, traj)

call getvalues_tl(geom, inc, locs, vars, gom, traj)

end subroutine fv3jedi_increment_getvalues_tl_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_getvalues_ad_c(c_key_geom, c_key_inc,c_key_loc,c_vars,c_key_gom,c_key_traj) bind(c,name='fv3jedi_increment_getvalues_ad_f90')
use iso_c_binding
use fv3jedi_increment_mod
use ioda_locs_mod
use ioda_locs_mod_c, only: ioda_locs_registry
use ufo_vars_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
use fv3jedi_geom_mod
use fv3jedi_getvaltraj_mod, only: fv3jedi_getvaltraj, fv3jedi_getvaltraj_registry
implicit none
integer(c_int), intent(in) :: c_key_inc  !< Increment to be interpolated
integer(c_int), intent(in) :: c_key_loc  !< List of requested locations
type(c_ptr), intent(in)    :: c_vars     !< List of requested variables
integer(c_int), intent(in) :: c_key_gom  !< Interpolated values
integer(c_int), intent(in) :: c_key_traj !< Trajectory for interpolation/transforms
integer(c_int), intent(in) :: c_key_geom  !< Geometry
type(fv3jedi_increment), pointer :: inc
type(ioda_locs),  pointer :: locs
type(ufo_geovals),  pointer :: gom
type(ufo_vars) :: vars
type(fv3jedi_getvaltraj), pointer :: traj
type(fv3jedi_geom),  pointer :: geom

call ufo_vars_setup(vars, c_vars)

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_increment_registry%get(c_key_inc, inc)
call ioda_locs_registry%get(c_key_loc, locs)
call ufo_geovals_registry%get(c_key_gom, gom)
call fv3jedi_getvaltraj_registry%get(c_key_traj, traj)

call getvalues_ad(geom, inc, locs, vars, gom, traj)

end subroutine fv3jedi_increment_getvalues_ad_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_sizes_c(c_key_self,nx,ny,nv) bind(c,name='fv3jedi_increment_sizes_f90')

use iso_c_binding
use fv3jedi_increment_mod

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: nx,ny,nv
type(fv3jedi_increment), pointer :: self

call fv3jedi_increment_registry%get(c_key_self,self)

nv = self%vars%nv
nx = self%npx
ny = self%npy

end subroutine fv3jedi_increment_sizes_c

! ------------------------------------------------------------------------------   
