! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_traj_interface_mod

use iso_c_binding
use fv3jedi_traj_mod
use fv3jedi_model_mod, only: fv3jedi_model
use fv3jedi_state_mod, only: fv3jedi_state, fv3jedi_state_registry
use fv3jedi_model_interface_mod, only: fv3jedi_model_registry

implicit none
private

public :: fv3jedi_traj_registry

!> Linked list interface - defines registry_t type
#define LISTED_TYPE fv3jedi_traj
#include "linkedList_i.f"
type(registry_t) :: fv3jedi_traj_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_traj_prop(c_key_model, c_key_state, c_key_self) bind(c,name='fv3jedi_traj_prop_f90')

implicit none
integer(c_int), intent(in)    :: c_key_model !< Model
integer(c_int), intent(in)    :: c_key_state !< State
integer(c_int), intent(inout) :: c_key_self  !< traj

type(fv3jedi_model), pointer :: model
type(fv3jedi_state), pointer :: state
type(fv3jedi_traj),  pointer :: self

call fv3jedi_model_registry%get(c_key_model,model)
call fv3jedi_state_registry%get(c_key_state,state)

call fv3jedi_traj_registry%init()            
call fv3jedi_traj_registry%add(c_key_self)
call fv3jedi_traj_registry%get(c_key_self,self)

call traj_prop(model, state, self)

end subroutine c_fv3jedi_traj_prop

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_traj_wipe(c_key_self) bind(c,name='fv3jedi_traj_wipe_f90')

implicit none
integer(c_int), intent(inout)   :: c_key_self  !< traj

type(fv3jedi_traj), pointer :: self

call fv3jedi_traj_registry%get(c_key_self,self)

call traj_wipe(self)

call fv3jedi_traj_registry%remove(c_key_self)

end subroutine c_fv3jedi_traj_wipe

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_traj_minmaxrms(c_key_self, pminmax) bind(c,name='fv3jedi_traj_minmaxrms_f90')

implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(inout) :: pminmax(3,11)

type(fv3jedi_traj), pointer :: self

call fv3jedi_traj_registry%get(c_key_self,self)

call traj_minmaxrms(self,pminmax)

end subroutine c_fv3jedi_traj_minmaxrms

! ------------------------------------------------------------------------------

end module fv3jedi_traj_interface_mod
