! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_model_setup(c_conf, c_key_geom, c_key_self) bind (c,name='fv3jedi_model_setup_f90')

use iso_c_binding
use config_mod
use duration_mod
use fv3jedi_model_mod
use fv3jedi_geom_mod
use kinds

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Key to model data
integer(c_int), intent(in)    :: c_key_geom  !< Geometry
type(c_ptr), intent(in)       :: c_conf      !< pointer to object of class Config

type(fv3jedi_model), pointer :: model
type(fv3jedi_geom), pointer :: geom

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_model_registry%init()
call fv3jedi_model_registry%add(c_key_self)
call fv3jedi_model_registry%get(c_key_self, model)

call model_setup(model, geom, c_conf)

end subroutine c_fv3jedi_model_setup

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_model_delete(c_key_self) bind (c,name='fv3jedi_model_delete_f90')

use fv3jedi_model_mod
use iso_c_binding

implicit none
integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_model), pointer :: self

call fv3jedi_model_registry%get(c_key_self, self)

call model_delete(self)

call fv3jedi_model_registry%remove(c_key_self)

end subroutine c_fv3jedi_model_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_model_prepare_integration(c_key_self, c_key_state) &
         & bind(c,name='fv3jedi_model_prepare_integration_f90')

use iso_c_binding
use fv3jedi_state_mod
use fv3jedi_model_mod

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state

type(fv3jedi_model), pointer :: self
type(fv3jedi_state), pointer :: state

call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_model_registry%get(c_key_self, self)

call model_prepare_integration(self, state)

end subroutine c_fv3jedi_model_prepare_integration

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_model_prepare_integration_ad(c_key_self, c_key_incr) &
           bind(c,name='fv3jedi_model_prepare_integration_ad_f90')

use iso_c_binding
use fv3jedi_increment_mod
use fv3jedi_model_mod

implicit none
integer(c_int), intent(in) :: c_key_self !< Model
integer(c_int), intent(in) :: c_key_incr !< Model increment

type(fv3jedi_model), pointer :: self
type(fv3jedi_increment), pointer :: inc

call fv3jedi_model_registry%get(c_key_self,self)
call fv3jedi_increment_registry%get(c_key_incr,inc)

call model_prepare_integration_ad(self, inc)

end subroutine c_fv3jedi_model_prepare_integration_ad

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_model_prepare_integration_tl(c_key_self, c_key_incr) &
           bind(c,name='fv3jedi_model_prepare_integration_tl_f90')

use iso_c_binding
use fv3jedi_increment_mod
use fv3jedi_model_mod

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_incr  !< Model increment

type(fv3jedi_model), pointer :: self
type(fv3jedi_increment), pointer :: inc

call fv3jedi_model_registry%get(c_key_self, self)
call fv3jedi_increment_registry%get(c_key_incr, inc)

call model_prepare_integration_tl(self, inc)

end subroutine c_fv3jedi_model_prepare_integration_tl

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_model_propagate(c_key_geom, c_key_self, c_key_state) bind(c,name='fv3jedi_model_propagate_f90')

use iso_c_binding
use fv3jedi_state_mod
use fv3jedi_model_mod
use fv3jedi_geom_mod

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state
integer(c_int), intent(in)    :: c_key_geom  !< Geometry

type(fv3jedi_geom), pointer :: geom
type(fv3jedi_model), pointer :: self
type(fv3jedi_state), pointer :: state

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_model_registry%get(c_key_self, self)
call fv3jedi_state_registry%get(c_key_state,state)

call model_propagate(geom, self, state)

end subroutine c_fv3jedi_model_propagate

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_model_propagate_ad(c_key_geom, c_key_self, c_key_incr, c_key_traj) &
           bind(c,name='fv3jedi_model_propagate_ad_f90')

use iso_c_binding
use fv3jedi_increment_mod
use fv3jedi_trajectories
use fv3jedi_model_mod
use fv3jedi_geom_mod

implicit none
integer(c_int), intent(in) :: c_key_self !< Model
integer(c_int), intent(in) :: c_key_incr !< Model increment
integer(c_int), intent(in) :: c_key_traj !< Trajectory structure
integer(c_int), intent(in)    :: c_key_geom  !< Geometry

type(fv3jedi_geom), pointer :: geom
type(fv3jedi_model),      pointer :: self
type(fv3jedi_increment),      pointer :: inc
type(fv3jedi_trajectory), pointer :: traj

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_model_registry%get(c_key_self,self)
call fv3jedi_increment_registry%get(c_key_incr,inc)
call fv3jedi_traj_registry%get(c_key_traj,traj)

call model_propagate_ad(geom, self, inc, traj)

end subroutine c_fv3jedi_model_propagate_ad

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_model_propagate_tl(c_key_geom, c_key_self, c_key_incr, c_key_traj) &
           bind(c,name='fv3jedi_model_propagate_tl_f90')

use iso_c_binding
use fv3jedi_increment_mod
use fv3jedi_trajectories
use fv3jedi_model_mod
use fv3jedi_geom_mod

implicit none
integer(c_int), intent(in) :: c_key_self !< Model
integer(c_int), intent(in) :: c_key_incr !< Model increment
integer(c_int), intent(in) :: c_key_traj !< Trajectory structure
integer(c_int), intent(in)    :: c_key_geom  !< Geometry

type(fv3jedi_geom), pointer :: geom
type(fv3jedi_model),      pointer :: self
type(fv3jedi_increment),      pointer :: inc
type(fv3jedi_trajectory), pointer :: traj

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_model_registry%get(c_key_self, self)
call fv3jedi_increment_registry%get(c_key_incr,inc)
call fv3jedi_traj_registry%get(c_key_traj,traj)

call model_propagate_tl(geom, self, inc, traj)

end subroutine c_fv3jedi_model_propagate_tl

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_model_prop_traj(c_key_self, c_key_state, c_key_traj) bind(c,name='fv3jedi_model_prop_traj_f90')

use iso_c_binding
use fv3jedi_state_mod
use fv3jedi_model_mod
use fv3jedi_trajectories

implicit none
integer(c_int), intent(in)    :: c_key_self  !< Model
integer(c_int), intent(in)    :: c_key_state !< Model state
integer(c_int), intent(inout) :: c_key_traj  !< Trajectory structure

type(fv3jedi_model),      pointer :: self
type(fv3jedi_state),      pointer :: state
type(fv3jedi_trajectory), pointer :: traj

call fv3jedi_model_registry%get(c_key_self,self)
call fv3jedi_state_registry%get(c_key_state,state)

call fv3jedi_traj_registry%init()            
call fv3jedi_traj_registry%add(c_key_traj)
call fv3jedi_traj_registry%get(c_key_traj,traj)

call model_prop_traj(self, state, traj)

end subroutine c_fv3jedi_model_prop_traj

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_model_wipe_traj(c_key_traj) bind(c,name='fv3jedi_model_wipe_traj_f90')

use iso_c_binding
use fv3jedi_model_mod
use fv3jedi_trajectories

implicit none
integer(c_int), intent(inout)   :: c_key_traj  !< Trajectory structure
type(fv3jedi_trajectory), pointer :: traj

call fv3jedi_traj_registry%get(c_key_traj,traj)

call model_wipe_traj(traj)

call fv3jedi_traj_registry%remove(c_key_traj)

end subroutine c_fv3jedi_model_wipe_traj

! ------------------------------------------------------------------------------
