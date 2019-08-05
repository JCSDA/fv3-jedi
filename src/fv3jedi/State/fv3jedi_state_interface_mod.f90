! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! ------------------------------------------------------------------------------

module fv3jedi_state_interface_mod

use fv3jedi_kinds_mod
use datetime_mod
use duration_mod
use iso_c_binding
use variables_mod
use fckit_configuration_module, only: fckit_configuration

use fv3jedi_state_mod
use fv3jedi_state_utils_mod, only: fv3jedi_state_registry
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_geom_interface_mod, only: fv3jedi_geom_registry
use fv3jedi_increment_utils_mod, only: fv3jedi_increment, fv3jedi_increment_registry

!GetValues
use ufo_locs_mod
use ufo_locs_mod_c, only: ufo_locs_registry
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry
use fv3jedi_getvalues_traj_mod, only: fv3jedi_getvalues_traj, fv3jedi_getvalues_traj_registry

private
public :: fv3jedi_state_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_create_c(c_key_self, c_key_geom, c_vars) bind(c,name='fv3jedi_state_create_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in) :: c_key_geom !< Geometry
type(c_ptr), intent(in)    :: c_vars     !< List of variables

type(fv3jedi_state), pointer :: self
type(fv3jedi_geom),  pointer :: geom
type(oops_vars)              :: vars
type(fckit_configuration)    :: f_conf

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_state_registry%init()
call fv3jedi_state_registry%add(c_key_self)
call fv3jedi_state_registry%get(c_key_self,self)

f_conf = fckit_configuration(c_vars)

call oops_vars_create(f_conf,vars)
call create(self, geom, vars)
call oops_vars_delete(vars)

end subroutine fv3jedi_state_create_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_delete_c(c_key_self) bind(c,name='fv3jedi_state_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_state), pointer :: self

call fv3jedi_state_registry%get(c_key_self,self)

call delete(self)

call fv3jedi_state_registry%remove(c_key_self)

end subroutine fv3jedi_state_delete_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_zero_c(c_key_self) bind(c,name='fv3jedi_state_zero_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
type(fv3jedi_state), pointer :: self

call fv3jedi_state_registry%get(c_key_self,self)
call zeros(self)

end subroutine fv3jedi_state_zero_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_copy_c(c_key_self,c_key_rhs) bind(c,name='fv3jedi_state_copy_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_state), pointer :: self
type(fv3jedi_state), pointer :: rhs
call fv3jedi_state_registry%get(c_key_self,self)
call fv3jedi_state_registry%get(c_key_rhs,rhs)

call copy(self, rhs)

end subroutine fv3jedi_state_copy_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_axpy_c(c_key_self,c_zz,c_key_rhs) bind(c,name='fv3jedi_state_axpy_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
integer(c_int), intent(in) :: c_key_rhs

type(fv3jedi_state), pointer :: self
type(fv3jedi_state), pointer :: rhs
real(kind=kind_real) :: zz

call fv3jedi_state_registry%get(c_key_self,self)
call fv3jedi_state_registry%get(c_key_rhs,rhs)
zz = c_zz

call axpy(self,zz,rhs)

end subroutine fv3jedi_state_axpy_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_add_incr_c(c_key_geom,c_key_self,c_key_rhs) bind(c,name='fv3jedi_state_add_incr_f90')

implicit none
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs
type(fv3jedi_geom), pointer :: geom
type(fv3jedi_state), pointer :: self
type(fv3jedi_increment), pointer :: rhs

call fv3jedi_state_registry%get(c_key_self,self)
call fv3jedi_increment_registry%get(c_key_rhs,rhs)
call fv3jedi_geom_registry%get(c_key_geom, geom)

call add_incr(geom,self,rhs)

end subroutine fv3jedi_state_add_incr_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_change_resol_c(c_key_state,c_key_geom,c_key_rhs,c_key_geom_rhs) bind(c,name='fv3jedi_state_change_resol_f90')

implicit none
integer(c_int), intent(in) :: c_key_state
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_rhs
integer(c_int), intent(in) :: c_key_geom_rhs

type(fv3jedi_state), pointer :: state, rhs
type(fv3jedi_geom),  pointer :: geom, geom_rhs

call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_state_registry%get(c_key_rhs,rhs)
call fv3jedi_geom_registry%get(c_key_geom_rhs, geom_rhs)

call change_resol(state,geom,rhs,geom_rhs)

end subroutine fv3jedi_state_change_resol_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_read_file_c(c_key_geom, c_key_state, c_conf, c_dt) bind(c,name='fv3jedi_state_read_file_f90')

implicit none
integer(c_int), intent(in) :: c_key_state  !< State
type(c_ptr), intent(in)    :: c_conf !< Configuration
type(c_ptr), intent(inout) :: c_dt   !< DateTime
integer(c_int), intent(in) :: c_key_geom  !< Geometry

type(fv3jedi_state), pointer :: state
type(datetime) :: fdate
type(fv3jedi_geom),  pointer :: geom

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_state_registry%get(c_key_state,state)
call c_f_datetime(c_dt, fdate)
call read_file(geom, state, c_conf, fdate)

end subroutine fv3jedi_state_read_file_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_analytic_init_c(c_key_state, c_key_geom, c_conf, c_dt) bind(c,name='fv3jedi_state_analytic_init_f90')

implicit none
integer(c_int), intent(in) :: c_key_state  !< State
integer(c_int), intent(in) :: c_key_geom  !< Geometry
type(c_ptr), intent(in)    :: c_conf !< Configuration
type(c_ptr), intent(inout) :: c_dt   !< DateTime

type(fv3jedi_state), pointer :: state
type(fv3jedi_geom), pointer :: geom
type(datetime) :: fdate

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_state_registry%get(c_key_state,state)
call c_f_datetime(c_dt, fdate)
call analytic_IC(state, geom, c_conf, fdate)

end subroutine fv3jedi_state_analytic_init_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_write_file_c(c_key_geom, c_key_state, c_conf, c_dt) bind(c,name='fv3jedi_state_write_file_f90')

implicit none
integer(c_int), intent(in) :: c_key_state  !< State
type(c_ptr), intent(in) :: c_conf !< Configuration
type(c_ptr), intent(in) :: c_dt   !< DateTime
integer(c_int), intent(in) :: c_key_geom  !< Geometry

type(fv3jedi_state), pointer :: state
type(datetime) :: fdate
type(fv3jedi_geom),  pointer :: geom

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_state_registry%get(c_key_state,state)
call c_f_datetime(c_dt, fdate)
call write_file(geom, state, c_conf, fdate)

end subroutine fv3jedi_state_write_file_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_gpnorm_c(c_key_state, kf, pstat) bind(c,name='fv3jedi_state_gpnorm_f90')

implicit none
integer(c_int), intent(in) :: c_key_state
integer(c_int), intent(in) :: kf
real(c_double), intent(inout) :: pstat(3*kf)

type(fv3jedi_state), pointer :: state
real(kind=kind_real) :: zstat(3, kf)
integer :: jj, js, jf

call fv3jedi_state_registry%get(c_key_state,state)

call gpnorm(state, kf, zstat)
jj=0
do jf = 1, kf
  do js = 1, 3
    jj=jj+1
    pstat(jj) = zstat(js,jf)
  enddo
enddo

end subroutine fv3jedi_state_gpnorm_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_print_c(c_key_self) bind(c,name='fv3jedi_state_print_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
type(fv3jedi_state), pointer :: self

call fv3jedi_state_registry%get(c_key_self,self)

call state_print(self)

end subroutine fv3jedi_state_print_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_rms_c(c_key_state, prms) bind(c,name='fv3jedi_state_rms_f90')

implicit none
integer(c_int), intent(in) :: c_key_state
real(c_double), intent(inout) :: prms

type(fv3jedi_state), pointer :: state
real(kind=kind_real) :: zz

call fv3jedi_state_registry%get(c_key_state,state)

call rms(state, zz)

prms = zz

end subroutine fv3jedi_state_rms_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_getvalues_notraj_c(c_key_geom, c_key_state,c_key_loc,c_vars,c_key_gom) bind(c,name='fv3jedi_state_getvalues_notraj_f90')

implicit none
integer(c_int), intent(in) :: c_key_state  !< State to be interpolated
integer(c_int), intent(in) :: c_key_loc  !< List of requested locations
type(c_ptr), intent(in)    :: c_vars     !< List of requested variables
integer(c_int), intent(in) :: c_key_gom  !< Interpolated values
integer(c_int), intent(in) :: c_key_geom  !< Geometry
type(fv3jedi_state), pointer :: state
type(ufo_locs),  pointer :: locs
type(ufo_geovals),  pointer :: gom
type(oops_vars) :: vars
type(fv3jedi_geom),  pointer :: geom

type(fckit_configuration)    :: f_conf

f_conf = fckit_configuration(c_vars)

call oops_vars_create(f_conf,vars)

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_state_registry%get(c_key_state, state)
call ufo_locs_registry%get(c_key_loc, locs)
call ufo_geovals_registry%get(c_key_gom, gom)

call getvalues(geom, state, locs, vars, gom)

call oops_vars_delete(vars)

end subroutine fv3jedi_state_getvalues_notraj_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_getvalues_c(c_key_geom, c_key_state,c_key_loc,c_vars,c_key_gom,c_key_traj) bind(c,name='fv3jedi_state_getvalues_f90')

implicit none
integer(c_int), intent(in) :: c_key_state  !< State to be interpolated
integer(c_int), intent(in) :: c_key_loc  !< List of requested locations
type(c_ptr), intent(in)    :: c_vars     !< List of requested variables
integer(c_int), intent(in) :: c_key_gom  !< Interpolated values
integer(c_int), intent(in) :: c_key_traj !< Trajectory for interpolation/transforms
integer(c_int), intent(in) :: c_key_geom  !< Geometry

type(fv3jedi_state), pointer :: state
type(ufo_locs),  pointer :: locs
type(ufo_geovals),  pointer :: gom
type(oops_vars) :: vars
type(fv3jedi_getvalues_traj), pointer :: traj
type(fv3jedi_geom),  pointer :: geom
type(fckit_configuration)    :: f_conf

f_conf = fckit_configuration(c_vars)

call oops_vars_create(f_conf,vars)

call fv3jedi_state_registry%get(c_key_state, state)
call fv3jedi_geom_registry%get(c_key_geom, geom)
call ufo_locs_registry%get(c_key_loc, locs)
call ufo_geovals_registry%get(c_key_gom, gom)
call fv3jedi_getvalues_traj_registry%get(c_key_traj, traj)

call getvalues(geom, state, locs, vars, gom, traj)

call oops_vars_delete(vars)

end subroutine fv3jedi_state_getvalues_c

! ------------------------------------------------------------------------------

subroutine fv3jedi_state_sizes_c(c_key_self,nx,ny,nf) bind(c,name='fv3jedi_state_sizes_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(inout) :: nx,ny,nf
type(fv3jedi_state), pointer :: self

call fv3jedi_state_registry%get(c_key_self,self)

nf = self%nf
nx = self%npx
ny = self%npy

end subroutine fv3jedi_state_sizes_c

! ------------------------------------------------------------------------------

end module fv3jedi_state_interface_mod
