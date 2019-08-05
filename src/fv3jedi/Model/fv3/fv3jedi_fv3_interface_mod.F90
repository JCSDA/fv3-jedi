! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_fv3_interface_mod

use fv3jedi_kinds_mod
use datetime_mod
use duration_mod
use iso_c_binding

use fv3jedi_fv3_mod
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_geom_interface_mod, only: fv3jedi_geom_registry
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry

implicit none
private

public :: fv3jedi_fv3_registry

!> Linked list interface
#define LISTED_TYPE fv3_model
#include "Utilities/linkedList_i.f"
type(registry_t) :: fv3jedi_fv3_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "Utilities/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_fv3_create(c_conf, c_key_geom, c_key_self) bind (c,name='fv3jedi_fv3_create_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Key to model data
integer(c_int), intent(in)    :: c_key_geom  !< Geometry
type(c_ptr),    intent(in)    :: c_conf      !< pointer to object of class Config

type(fv3_model), pointer :: self
type(fv3jedi_geom),  pointer :: geom

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_fv3_registry%init()
call fv3jedi_fv3_registry%add(c_key_self)
call fv3jedi_fv3_registry%get(c_key_self, self)

call fv3_create(self, geom, c_conf)

end subroutine c_fv3jedi_fv3_create

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_fv3_delete(c_key_self) bind (c,name='fv3jedi_fv3_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(fv3_model), pointer :: self

call fv3jedi_fv3_registry%get(c_key_self, self)

call fv3_delete(self)

call fv3jedi_fv3_registry%remove(c_key_self)

end subroutine c_fv3jedi_fv3_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_fv3_initialize(c_key_self, c_key_state) bind(c,name='fv3jedi_fv3_initialize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state

type(fv3_model), pointer :: self
type(fv3jedi_state), pointer :: state

call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_fv3_registry%get(c_key_self, self)

call fv3_initialize(self, state)

end subroutine c_fv3jedi_fv3_initialize

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_fv3_step(c_key_self, c_key_state, c_key_geom, c_dt) bind(c,name='fv3jedi_fv3_step_f90')

implicit none
integer(c_int), intent(in)    :: c_key_self  !< Model
integer(c_int), intent(in)    :: c_key_state !< Model state
integer(c_int), intent(in)    :: c_key_geom  !< Geometry
type(c_ptr),    intent(inout) :: c_dt        !< DateTime

type(fv3_model),     pointer :: self
type(fv3jedi_state), pointer :: state
type(fv3jedi_geom),  pointer :: geom
type(datetime)               :: fdate

call fv3jedi_fv3_registry%get(c_key_self, self)
call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_geom_registry%get(c_key_geom, geom)
call c_f_datetime(c_dt, fdate)

call fv3_step(self, state, geom, fdate)

end subroutine c_fv3jedi_fv3_step

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_fv3_finalize(c_key_self, c_key_state) bind(c,name='fv3jedi_fv3_finalize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state

type(fv3_model), pointer :: self
type(fv3jedi_state), pointer :: state

call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_fv3_registry%get(c_key_self, self)

call fv3_finalize(self, state)

end subroutine c_fv3jedi_fv3_finalize

! ------------------------------------------------------------------------------

end module fv3jedi_fv3_interface_mod
