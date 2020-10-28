! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_traj_interface_mod

! iso
use iso_c_binding

! fv3-jedi
use fv3jedi_traj_mod,            only: fv3jedi_traj, wipe, set
use fv3jedi_state_mod,           only: fv3jedi_state
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry

implicit none
private

public :: fv3jedi_traj_registry

!> Linked list interface - defines registry_t type
#define LISTED_TYPE fv3jedi_traj
#include "oops/util/linkedList_i.f"
type(registry_t) :: fv3jedi_traj_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_traj_set(c_key_self, c_key_state) bind(c,name='fv3jedi_traj_set_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self  !< traj
integer(c_int), intent(in)    :: c_key_state !< State

type(fv3jedi_state), pointer :: state
type(fv3jedi_traj),  pointer :: self

! LinkedList
call fv3jedi_state_registry%get(c_key_state, state)
call fv3jedi_traj_registry%init()
call fv3jedi_traj_registry%add(c_key_self)
call fv3jedi_traj_registry%get(c_key_self, self)

! Implementation
call set(self, state)

end subroutine c_fv3jedi_traj_set

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_traj_wipe(c_key_self) bind(c,name='fv3jedi_traj_wipe_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(fv3jedi_traj), pointer :: self

! LinkedList
call fv3jedi_traj_registry%get(c_key_self,self)

! Implementation
call wipe(self)

! LinkedList
call fv3jedi_traj_registry%remove(c_key_self)

end subroutine c_fv3jedi_traj_wipe

! --------------------------------------------------------------------------------------------------

end module fv3jedi_traj_interface_mod
