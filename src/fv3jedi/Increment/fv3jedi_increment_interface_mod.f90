! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! ------------------------------------------------------------------------------

module fv3jedi_increment_interface_mod

use fv3jedi_kinds_mod
use datetime_mod
use duration_mod
use iso_c_binding
use variables_mod
use fckit_configuration_module, only: fckit_configuration

use fv3jedi_increment_mod
use fv3jedi_increment_utils_mod, only: fv3jedi_increment_registry
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_geom_interface_mod, only: fv3jedi_geom_registry
use fv3jedi_state_utils_mod, only: fv3jedi_state, fv3jedi_state_registry
use unstructured_grid_mod, only: unstructured_grid, unstructured_grid_registry

!GetValues
use ufo_locs_mod
use ufo_locs_mod_c, only: ufo_locs_registry
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
use fv3jedi_getvalues_traj_mod, only: fv3jedi_getvalues_traj, fv3jedi_getvalues_traj_registry

implicit none
private

public :: fv3jedi_increment_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_create_c(c_key_self, c_key_geom, c_vars) bind(c,name='fv3jedi_increment_create_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in) :: c_key_geom !< Geometry
type(c_ptr), intent(in)    :: c_vars     !< List of variables

type(fv3jedi_increment), pointer :: self
type(fv3jedi_geom),  pointer :: geom
type(oops_vars) :: vars
type(fckit_configuration)    :: f_conf

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_increment_registry%init()
call fv3jedi_increment_registry%add(c_key_self)
call fv3jedi_increment_registry%get(c_key_self,self)

f_conf = fckit_configuration(c_vars)

call oops_vars_create(f_conf,vars)
call create(self, geom, vars)
call oops_vars_delete(vars)

end subroutine fv3jedi_increment_create_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_delete_c(c_key_self) bind(c,name='fv3jedi_increment_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_increment), pointer :: self

call fv3jedi_increment_registry%get(c_key_self,self)

call delete(self)

call fv3jedi_increment_registry%remove(c_key_self)

end subroutine fv3jedi_increment_delete_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_zero_c(c_key_self) bind(c,name='fv3jedi_increment_zero_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
type(fv3jedi_increment), pointer :: self

call fv3jedi_increment_registry%get(c_key_self,self)
call zeros(self)

end subroutine fv3jedi_increment_zero_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_dirac_c(c_key_self,c_conf,c_key_geom) bind(c,name='fv3jedi_increment_dirac_f90')

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

implicit none
integer(c_int), intent(in) :: c_key_self
type(fv3jedi_increment), pointer :: self

call fv3jedi_increment_registry%get(c_key_self,self)
call random(self)

end subroutine fv3jedi_increment_random_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_ug_coord_c(c_key_inc, c_key_ug, c_key_geom) bind (c,name='fv3jedi_increment_ug_coord_f90')

implicit none
integer(c_int), intent(in) :: c_key_inc
integer(c_int), intent(in) :: c_key_ug
integer(c_int), intent(in) :: c_key_geom !< Geometry
type(fv3jedi_increment), pointer :: inc
type(unstructured_grid), pointer :: ug
type(fv3jedi_geom),  pointer :: geom

call fv3jedi_increment_registry%get(c_key_inc,inc)
call unstructured_grid_registry%get(c_key_ug,ug)
call fv3jedi_geom_registry%get(c_key_geom, geom)

call ug_coord(inc, ug, geom)

end subroutine fv3jedi_increment_ug_coord_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_increment_to_ug_c(c_key_inc, c_key_ug, c_its) bind (c,name='fv3jedi_increment_increment_to_ug_f90')

implicit none
integer(c_int), intent(in) :: c_key_inc
integer(c_int), intent(in) :: c_key_ug
integer(c_int), intent(in) :: c_its
type(fv3jedi_increment), pointer :: inc
type(unstructured_grid), pointer :: ug
integer :: its

its = c_its+1

call fv3jedi_increment_registry%get(c_key_inc,inc)
call unstructured_grid_registry%get(c_key_ug,ug)

call increment_to_ug(inc, ug, its)

end subroutine fv3jedi_increment_increment_to_ug_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_increment_from_ug_c(c_key_inc, c_key_ug, c_its) bind (c,name='fv3jedi_increment_increment_from_ug_f90')

implicit none
integer(c_int), intent(in) :: c_key_inc
integer(c_int), intent(in) :: c_key_ug
integer(c_int), intent(in) :: c_its
type(fv3jedi_increment), pointer :: inc
type(unstructured_grid), pointer :: ug
integer :: its

its = c_its+1

call fv3jedi_increment_registry%get(c_key_inc,inc)
call unstructured_grid_registry%get(c_key_ug,ug)

call increment_from_ug(inc, ug, its)

end subroutine fv3jedi_increment_increment_from_ug_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_copy_c(c_key_self,c_key_rhs) bind(c,name='fv3jedi_increment_copy_f90')

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

subroutine fv3jedi_increment_diff_incr_c(c_key_lhs,c_key_x1,c_key_x2,c_key_geom) bind(c,name='fv3jedi_increment_diff_incr_f90')

implicit none
integer(c_int), intent(in) :: c_key_lhs
integer(c_int), intent(in) :: c_key_x1
integer(c_int), intent(in) :: c_key_x2
integer(c_int), intent(in) :: c_key_geom

type(fv3jedi_increment), pointer :: lhs
type(fv3jedi_state), pointer :: x1
type(fv3jedi_state), pointer :: x2
type(fv3jedi_geom),  pointer :: geom

call fv3jedi_increment_registry%get(c_key_lhs,lhs)
call fv3jedi_state_registry%get(c_key_x1,x1)
call fv3jedi_state_registry%get(c_key_x2,x2)
call fv3jedi_geom_registry%get(c_key_geom, geom)

call diff_incr(lhs,x1,x2,geom)

end subroutine fv3jedi_increment_diff_incr_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_change_resol_c(c_key_inc,c_key_geom,c_key_rhs,c_key_geom_rhs) bind(c,name='fv3jedi_increment_change_resol_f90')

implicit none
integer(c_int), intent(in) :: c_key_inc
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_rhs
integer(c_int), intent(in) :: c_key_geom_rhs

type(fv3jedi_increment), pointer :: inc, rhs
type(fv3jedi_geom),  pointer :: geom, geom_rhs

call fv3jedi_increment_registry%get(c_key_inc,inc)
call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_increment_registry%get(c_key_rhs,rhs)
call fv3jedi_geom_registry%get(c_key_geom_rhs, geom_rhs)

call change_resol(inc,geom,rhs,geom_rhs)

end subroutine fv3jedi_increment_change_resol_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_read_file_c(c_key_geom, c_key_inc, c_conf, c_dt) bind(c,name='fv3jedi_increment_read_file_f90')

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

subroutine fv3jedi_increment_print_c(c_key_self) bind(c,name='fv3jedi_increment_print_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
type(fv3jedi_increment), pointer :: self

call fv3jedi_increment_registry%get(c_key_self,self)

call increment_print(self)

end subroutine fv3jedi_increment_print_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_rms_c(c_key_inc, prms) bind(c,name='fv3jedi_increment_rms_f90')

implicit none
integer(c_int), intent(in) :: c_key_inc
real(c_double), intent(inout) :: prms

type(fv3jedi_increment), pointer :: inc
real(kind=kind_real) :: zz

call fv3jedi_increment_registry%get(c_key_inc,inc)

call rms(inc, zz)

prms = zz

end subroutine fv3jedi_increment_rms_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_getvalues_tl_c(c_key_geom, c_key_inc,c_key_loc,c_vars,c_key_gom,c_key_traj) bind(c,name='fv3jedi_increment_getvalues_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_inc  !< Increment to be interpolated
integer(c_int), intent(in) :: c_key_loc  !< List of requested locations
type(c_ptr), intent(in)    :: c_vars     !< List of requested variables
integer(c_int), intent(in) :: c_key_gom  !< Interpolated values
integer(c_int), intent(in) :: c_key_traj !< Trajectory for interpolation/transforms
integer(c_int), intent(in) :: c_key_geom  !< Geometry
type(fv3jedi_increment), pointer :: inc
type(ufo_locs),  pointer :: locs
type(ufo_geovals),  pointer :: gom
type(oops_vars) :: vars
type(fv3jedi_getvalues_traj), pointer :: traj
type(fv3jedi_geom),  pointer :: geom
type(fckit_configuration)    :: f_conf

f_conf = fckit_configuration(c_vars)

call oops_vars_create(f_conf,vars)

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_increment_registry%get(c_key_inc, inc)
call ufo_locs_registry%get(c_key_loc, locs)
call ufo_geovals_registry%get(c_key_gom, gom)
call fv3jedi_getvalues_traj_registry%get(c_key_traj, traj)

call getvalues_tl(geom, inc, locs, vars, gom, traj)

call oops_vars_delete(vars)

end subroutine fv3jedi_increment_getvalues_tl_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_getvalues_ad_c(c_key_geom, c_key_inc,c_key_loc,c_vars,c_key_gom,c_key_traj) bind(c,name='fv3jedi_increment_getvalues_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_inc  !< Increment to be interpolated
integer(c_int), intent(in) :: c_key_loc  !< List of requested locations
type(c_ptr), intent(in)    :: c_vars     !< List of requested variables
integer(c_int), intent(in) :: c_key_gom  !< Interpolated values
integer(c_int), intent(in) :: c_key_traj !< Trajectory for interpolation/transforms
integer(c_int), intent(in) :: c_key_geom  !< Geometry
type(fv3jedi_increment), pointer :: inc
type(ufo_locs),  pointer :: locs
type(ufo_geovals),  pointer :: gom
type(oops_vars) :: vars
type(fv3jedi_getvalues_traj), pointer :: traj
type(fv3jedi_geom),  pointer :: geom
type(fckit_configuration)    :: f_conf

f_conf = fckit_configuration(c_vars)

call oops_vars_create(f_conf,vars)

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_increment_registry%get(c_key_inc, inc)
call ufo_locs_registry%get(c_key_loc, locs)
call ufo_geovals_registry%get(c_key_gom, gom)
call fv3jedi_getvalues_traj_registry%get(c_key_traj, traj)

call getvalues_ad(geom, inc, locs, vars, gom, traj)

call oops_vars_delete(vars)

end subroutine fv3jedi_increment_getvalues_ad_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_sizes_c(c_key_self,nx,ny,nv) bind(c,name='fv3jedi_increment_sizes_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: nx,ny,nv
type(fv3jedi_increment), pointer :: self

call fv3jedi_increment_registry%get(c_key_self,self)

nv = self%nf
nx = self%npx
ny = self%npy

end subroutine fv3jedi_increment_sizes_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_increment_jnormgrad_c(c_key_self,c_key_geom,c_key_state,c_conf) bind(c,name='fv3jedi_increment_jnormgrad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_state
type(c_ptr),    intent(in) :: c_conf

type(fv3jedi_increment), pointer :: self
type(fv3jedi_geom), pointer :: geom
type(fv3jedi_state), pointer :: state

call fv3jedi_increment_registry%get(c_key_self,self)
call fv3jedi_geom_registry%get(c_key_geom,geom)
call fv3jedi_state_registry%get(c_key_state,state)

call jnormgrad(self,geom,state,c_conf)

end subroutine fv3jedi_increment_jnormgrad_c

! ------------------------------------------------------------------------------
end module fv3jedi_increment_interface_mod
