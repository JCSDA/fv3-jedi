! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_fv3_interface_mod

use iso_c_binding

! oops uses
use datetime_mod
use duration_mod

! fckit uses
use fckit_configuration_module, only: fckit_configuration

! fv3-jedi uses
use fv3jedi_fv3_mod,             only: fv3_model
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_geom_interface_mod,  only: fv3jedi_geom_registry
use fv3jedi_state_mod,           only: fv3jedi_state
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry

implicit none
private
public :: fv3jedi_fv3_registry

!> Linked list interface
#define LISTED_TYPE fv3_model
#include "oops/util/linkedList_i.f"
type(registry_t) :: fv3jedi_fv3_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_fv3_create(c_conf, c_key_geom, c_key_self) &
           bind (c,name='fv3jedi_fv3_create_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Key to model data
integer(c_int), intent(in)    :: c_key_geom  !< Geometry
type(c_ptr),    intent(in)    :: c_conf      !< pointer to object of class Config

type(fv3_model),    pointer :: self
type(fv3jedi_geom), pointer :: geom
type(fckit_configuration)   :: f_conf

! Linked lists
call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_fv3_registry%init()
call fv3jedi_fv3_registry%add(c_key_self)
call fv3jedi_fv3_registry%get(c_key_self, self)

! Fortran configuration
f_conf = fckit_configuration(c_conf)

! Implementation
call self%create(geom, f_conf)

end subroutine c_fv3jedi_fv3_create

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_fv3_delete(c_key_self) bind (c,name='fv3jedi_fv3_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(fv3_model), pointer :: self

! Linked lists
call fv3jedi_fv3_registry%get(c_key_self, self)

! Implementation
call self%delete()

! Linked lists
call fv3jedi_fv3_registry%remove(c_key_self)

end subroutine c_fv3jedi_fv3_delete

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_fv3_initialize(c_key_self, c_key_state, c_dt) &
           bind(c,name='fv3jedi_fv3_initialize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state
type(c_ptr),    intent(inout) :: c_dt        !< DateTime

type(fv3_model),     pointer :: self
type(fv3jedi_state), pointer :: state
type(datetime)               :: fdate

! Linked lists
call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_fv3_registry%get(c_key_self, self)
call c_f_datetime(c_dt, fdate)

! Implementation
call self%initialize(state, fdate)

end subroutine c_fv3jedi_fv3_initialize

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_fv3_step(c_key_self, c_key_state, c_dt) &
           bind(c,name='fv3jedi_fv3_step_f90')

implicit none
integer(c_int), intent(in)    :: c_key_self  !< Model
integer(c_int), intent(in)    :: c_key_state !< Model state
type(c_ptr),    intent(inout) :: c_dt        !< DateTime

type(fv3_model),     pointer :: self
type(fv3jedi_state), pointer :: state
type(datetime)               :: fdate

! Linked lists
call fv3jedi_fv3_registry%get(c_key_self, self)
call fv3jedi_state_registry%get(c_key_state,state)
call c_f_datetime(c_dt, fdate)

! Implementation
call self%step(state, fdate)

end subroutine c_fv3jedi_fv3_step

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_fv3_finalize(c_key_self, c_key_state) bind(c,name='fv3jedi_fv3_finalize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state

type(fv3_model),     pointer :: self
type(fv3jedi_state), pointer :: state

! Linked lists
call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_fv3_registry%get(c_key_self, self)

! Implementation
call self%finalize(state)

end subroutine c_fv3jedi_fv3_finalize

! --------------------------------------------------------------------------------------------------

end module fv3jedi_fv3_interface_mod
