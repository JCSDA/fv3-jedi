! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_geos_interface_mod

use fv3jedi_kinds_mod
use datetime_mod
use duration_mod
use iso_c_binding

use fckit_configuration_module, only: fckit_configuration

use fv3jedi_geos_mod
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_geom_interface_mod, only: fv3jedi_geom_registry
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry

implicit none
private

public :: fv3jedi_geos_registry

!> Linked list interface
#define LISTED_TYPE geos_model
#include "oops/util/linkedList_i.f"
type(registry_t) :: fv3jedi_geos_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geos_create(c_conf, c_key_geom, c_key_self) bind (c,name='fv3jedi_geos_create_f90')

implicit none
integer(c_int), intent(inout)  :: c_key_self  !< Key to model data
integer(c_int), intent(in)     :: c_key_geom  !< Geometry
type(c_ptr), value, intent(in) :: c_conf      !< pointer to object of class Config

type(geos_model), pointer   :: self
type(fv3jedi_geom), pointer :: geom
type(fckit_configuration)   :: f_conf

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_geos_registry%init()
call fv3jedi_geos_registry%add(c_key_self)
call fv3jedi_geos_registry%get(c_key_self, self)

f_conf = fckit_configuration(c_conf)

call geos_create(self, geom, f_conf)

end subroutine c_fv3jedi_geos_create

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geos_delete(c_key_self) bind (c,name='fv3jedi_geos_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(geos_model), pointer :: self

call fv3jedi_geos_registry%get(c_key_self, self)

call geos_delete(self)

call fv3jedi_geos_registry%remove(c_key_self)

end subroutine c_fv3jedi_geos_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geos_initialize(c_key_self, c_key_state) bind(c,name='fv3jedi_geos_initialize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state

type(geos_model), pointer :: self
type(fv3jedi_state), pointer :: state

call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_geos_registry%get(c_key_self, self)

call geos_initialize(self, state)

end subroutine c_fv3jedi_geos_initialize

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geos_step(c_key_self, c_key_state) bind(c,name='fv3jedi_geos_step_f90')

implicit none
integer(c_int), intent(in)    :: c_key_self  !< Model
integer(c_int), intent(in)    :: c_key_state !< Model state

type(geos_model), pointer :: self
type(fv3jedi_state), pointer :: state

call fv3jedi_geos_registry%get(c_key_self, self)
call fv3jedi_state_registry%get(c_key_state,state)

call geos_step(self, state)

end subroutine c_fv3jedi_geos_step

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geos_finalize(c_key_self, c_key_state) bind(c,name='fv3jedi_geos_finalize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state

type(geos_model), pointer :: self
type(fv3jedi_state), pointer :: state

call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_geos_registry%get(c_key_self, self)

call geos_finalize(self, state)

end subroutine c_fv3jedi_geos_finalize

! ------------------------------------------------------------------------------

end module fv3jedi_geos_interface_mod
