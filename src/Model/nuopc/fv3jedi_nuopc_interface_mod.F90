! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_nuopc_interface_mod

use fv3jedi_kinds_mod
use config_mod
use datetime_mod
use duration_mod
use iso_c_binding

use fv3jedi_nuopc_mod
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_geom_interface_mod, only: fv3jedi_geom_registry
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry

implicit none
private

public :: fv3jedi_nuopc_registry

!> Linked list interface
#define LISTED_TYPE model_nuopc_type
#include "linkedList_i.f"
type(registry_t) :: fv3jedi_nuopc_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_nuopc_create(c_conf, c_key_geom, c_key_self) bind (c,name='fv3jedi_nuopc_create_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Key to model data
integer(c_int), intent(in)    :: c_key_geom  !< Geometry
type(c_ptr),    intent(in)    :: c_conf      !< pointer to object of class Config

type(model_nuopc_type), pointer :: self
type(fv3jedi_geom),  pointer :: geom

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_nuopc_registry%init()
call fv3jedi_nuopc_registry%add(c_key_self)
call fv3jedi_nuopc_registry%get(c_key_self, self)

call model_nuopc_create(self, geom, c_conf)

end subroutine c_fv3jedi_nuopc_create

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_nuopc_delete(c_key_self) bind (c,name='fv3jedi_nuopc_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(model_nuopc_type), pointer :: self

call fv3jedi_nuopc_registry%get(c_key_self, self)

call model_nuopc_delete(self)

call fv3jedi_nuopc_registry%remove(c_key_self)

end subroutine c_fv3jedi_nuopc_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_nuopc_initialize(c_key_self, c_key_state, c_dt) bind(c,name='fv3jedi_nuopc_initialize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state

type(model_nuopc_type), pointer :: self
type(fv3jedi_state), pointer :: state
type(c_ptr),    intent(inout) :: c_dt        !< DateTime

type(datetime) :: fdate

call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_nuopc_registry%get(c_key_self, self)

call c_f_datetime(c_dt, fdate)

call model_nuopc_initialize(self, state, fdate)

end subroutine c_fv3jedi_nuopc_initialize

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_nuopc_step(c_key_self, c_key_state, c_dt1, c_dt2) bind(c,name='fv3jedi_nuopc_step_f90')

implicit none
integer(c_int), intent(in)    :: c_key_self  !< Model
integer(c_int), intent(in)    :: c_key_state !< Model state
type(c_ptr),    intent(inout) :: c_dt1       !< DateTime
type(c_ptr),    intent(inout) :: c_dt2       !< DateTime

type(model_nuopc_type), pointer :: self
type(fv3jedi_state), pointer :: state

type(datetime) :: fdate1
type(datetime) :: fdate2

call fv3jedi_nuopc_registry%get(c_key_self, self)
call fv3jedi_state_registry%get(c_key_state,state)

call c_f_datetime(c_dt1, fdate1)
call c_f_datetime(c_dt2, fdate2)

call model_nuopc_step(self, state, fdate1, fdate2)

end subroutine c_fv3jedi_nuopc_step

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_nuopc_finalize(c_key_self, c_key_state, c_dt) bind(c,name='fv3jedi_nuopc_finalize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state
type(c_ptr),    intent(inout) :: c_dt        !< DateTime

type(model_nuopc_type), pointer :: self
type(fv3jedi_state), pointer :: state

type(datetime) :: fdate

call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_nuopc_registry%get(c_key_self, self)

call c_f_datetime(c_dt, fdate)

call model_nuopc_finalize(self, state, fdate)

end subroutine c_fv3jedi_nuopc_finalize

! ------------------------------------------------------------------------------

end module fv3jedi_nuopc_interface_mod
