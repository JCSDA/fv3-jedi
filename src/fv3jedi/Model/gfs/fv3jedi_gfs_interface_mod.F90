! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_gfs_interface_mod

use fv3jedi_kinds_mod
use config_mod
use datetime_mod
use duration_mod
use iso_c_binding

use fv3jedi_gfs_mod
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_geom_interface_mod, only: fv3jedi_geom_registry
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry

implicit none
private

public :: fv3jedi_gfs_registry

!> Linked list interface
#define LISTED_TYPE gfs_model
#include "Utilities/linkedList_i.f"
type(registry_t) :: fv3jedi_gfs_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "Utilities/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_gfs_create(c_conf, c_key_geom, c_key_self) bind (c,name='fv3jedi_gfs_create_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Key to model data
integer(c_int), intent(in)    :: c_key_geom  !< Geometry
type(c_ptr),    intent(in)    :: c_conf      !< pointer to object of class Config

type(gfs_model), pointer :: self
type(fv3jedi_geom),  pointer :: geom

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_gfs_registry%init()
call fv3jedi_gfs_registry%add(c_key_self)
call fv3jedi_gfs_registry%get(c_key_self, self)

call gfs_create(self, geom, c_conf)

end subroutine c_fv3jedi_gfs_create

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_gfs_delete(c_key_self) bind (c,name='fv3jedi_gfs_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(gfs_model), pointer :: self

call fv3jedi_gfs_registry%get(c_key_self, self)

call gfs_delete(self)

call fv3jedi_gfs_registry%remove(c_key_self)

end subroutine c_fv3jedi_gfs_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_gfs_initialize(c_key_self, c_key_state) bind(c,name='fv3jedi_gfs_initialize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state

type(gfs_model), pointer :: self
type(fv3jedi_state), pointer :: state

call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_gfs_registry%get(c_key_self, self)

call gfs_initialize(self, state)

end subroutine c_fv3jedi_gfs_initialize

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_gfs_step(c_key_self, c_key_state, c_dt) bind(c,name='fv3jedi_gfs_step_f90')

implicit none
integer(c_int), intent(in)    :: c_key_self  !< Model
integer(c_int), intent(in)    :: c_key_state !< Model state
type(c_ptr),    intent(inout) :: c_dt        !< DateTime

type(gfs_model), pointer :: self
type(fv3jedi_state), pointer :: state

type(datetime) :: fdate

call fv3jedi_gfs_registry%get(c_key_self, self)
call fv3jedi_state_registry%get(c_key_state,state)

call c_f_datetime(c_dt, fdate)

call gfs_step(self, state, fdate)

end subroutine c_fv3jedi_gfs_step

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_gfs_finalize(c_key_self, c_key_state) bind(c,name='fv3jedi_gfs_finalize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state

type(gfs_model), pointer :: self
type(fv3jedi_state), pointer :: state

call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_gfs_registry%get(c_key_self, self)

call gfs_finalize(self, state)

end subroutine c_fv3jedi_gfs_finalize

! ------------------------------------------------------------------------------

end module fv3jedi_gfs_interface_mod
