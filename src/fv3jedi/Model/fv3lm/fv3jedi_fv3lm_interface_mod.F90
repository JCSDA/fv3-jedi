! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_fv3lm_interface_mod

use iso_c_binding

! oops uses
use datetime_mod
use duration_mod

! fckit uses
use fckit_configuration_module, only: fckit_configuration

! fv3-jedi uses
use fv3jedi_fv3lm_mod,           only: fv3lm_model
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_geom_interface_mod,  only: fv3jedi_geom_registry
use fv3jedi_state_mod,           only: fv3jedi_state
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry

implicit none
private
public :: fv3jedi_fv3lm_registry

!> Linked list interface
#define LISTED_TYPE fv3lm_model
#include "oops/util/linkedList_i.f"
type(registry_t) :: fv3jedi_fv3lm_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_fv3lm_create(c_conf, c_key_geom, c_key_self) &
           bind (c,name='fv3jedi_fv3lm_create_f90')

implicit none
integer(c_int), intent(inout)  :: c_key_self  !< Key to model data
integer(c_int), intent(in)     :: c_key_geom  !< Geometry
type(c_ptr),value , intent(in) :: c_conf      !< pointer to object of class Config

type(fv3lm_model),  pointer :: self
type(fv3jedi_geom), pointer :: geom
type(fckit_configuration)   :: f_conf

! Linked lists
call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_fv3lm_registry%init()
call fv3jedi_fv3lm_registry%add(c_key_self)
call fv3jedi_fv3lm_registry%get(c_key_self, self)

! Fortran configuration
f_conf = fckit_configuration(c_conf)

! Implementation
call self%create(geom, f_conf)

end subroutine c_fv3jedi_fv3lm_create

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_fv3lm_delete(c_key_self) bind (c,name='fv3jedi_fv3lm_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(fv3lm_model), pointer :: self

! Linked lists
call fv3jedi_fv3lm_registry%get(c_key_self, self)

! Implementation
call self%delete()

! Linked lists
call fv3jedi_fv3lm_registry%remove(c_key_self)

end subroutine c_fv3jedi_fv3lm_delete

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_fv3lm_initialize(c_key_self, c_key_state) &
           bind(c,name='fv3jedi_fv3lm_initialize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state

type(fv3lm_model),   pointer :: self
type(fv3jedi_state), pointer :: state

! Linked lists
call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_fv3lm_registry%get(c_key_self, self)

! Implementation
call self%initialize(state)

end subroutine c_fv3jedi_fv3lm_initialize

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_fv3lm_step(c_key_self, c_key_state, c_key_geom) &
           bind(c,name='fv3jedi_fv3lm_step_f90')

implicit none
integer(c_int), intent(in)    :: c_key_self  !< Model
integer(c_int), intent(in)    :: c_key_state !< Model state
integer(c_int), intent(in)    :: c_key_geom  !< Geometry

type(fv3lm_model),   pointer :: self
type(fv3jedi_state), pointer :: state
type(fv3jedi_geom),  pointer :: geom

! Linked lists
call fv3jedi_fv3lm_registry%get(c_key_self, self)
call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_geom_registry%get(c_key_geom, geom)

! Implementation
call self%step(state, geom)

end subroutine c_fv3jedi_fv3lm_step

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_fv3lm_finalize(c_key_self, c_key_state) &
          bind(c,name='fv3jedi_fv3lm_finalize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state

type(fv3lm_model),   pointer :: self
type(fv3jedi_state), pointer :: state

! Linked lists
call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_fv3lm_registry%get(c_key_self, self)

! Implementation
call self%finalize(state)

end subroutine c_fv3jedi_fv3lm_finalize

! --------------------------------------------------------------------------------------------------

end module fv3jedi_fv3lm_interface_mod
