! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_model_interface_mod

use kinds
use config_mod
use duration_mod
use iso_c_binding

use fv3jedi_model_mod
use fv3jedi_geom_mod, only: fv3jedi_geom, fv3jedi_geom_registry
use fv3jedi_state_mod, only: fv3jedi_state, fv3jedi_state_registry

implicit none
private

public :: fv3jedi_model_registry

!> Linked list interface
#define LISTED_TYPE fv3jedi_model
#include "linkedList_i.f"
type(registry_t) :: fv3jedi_model_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_model_create(c_conf, c_key_geom, c_key_self) bind (c,name='fv3jedi_model_create_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Key to model data
integer(c_int), intent(in)    :: c_key_geom  !< Geometry
type(c_ptr),    intent(in)    :: c_conf      !< pointer to object of class Config

type(fv3jedi_model), pointer :: self
type(fv3jedi_geom),  pointer :: geom

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_model_registry%init()
call fv3jedi_model_registry%add(c_key_self)
call fv3jedi_model_registry%get(c_key_self, self)

call model_create(self, geom, c_conf)

end subroutine c_fv3jedi_model_create

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_model_delete(c_key_self) bind (c,name='fv3jedi_model_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_model), pointer :: self

call fv3jedi_model_registry%get(c_key_self, self)

call model_delete(self)

call fv3jedi_model_registry%remove(c_key_self)

end subroutine c_fv3jedi_model_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_model_initialize(c_key_self, c_key_state) bind(c,name='fv3jedi_model_initialize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state

type(fv3jedi_model), pointer :: self
type(fv3jedi_state), pointer :: state

call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_model_registry%get(c_key_self, self)

call model_initialize(self, state)

end subroutine c_fv3jedi_model_initialize

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_model_step(c_key_geom, c_key_self, c_key_state) bind(c,name='fv3jedi_model_step_f90')

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

call model_step(geom, self, state)

end subroutine c_fv3jedi_model_step

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_model_finalize(c_key_self, c_key_state) bind(c,name='fv3jedi_model_finalize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state

type(fv3jedi_model), pointer :: self
type(fv3jedi_state), pointer :: state

call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_model_registry%get(c_key_self, self)

call model_finalize(self, state)

end subroutine c_fv3jedi_model_finalize

! ------------------------------------------------------------------------------

end module fv3jedi_model_interface_mod
