! (C) Copyright 2017-2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module fv3jedi_tlm_interface_mod

use kinds
use config_mod
use duration_mod
use iso_c_binding

use fv3jedi_tlm_mod
use fv3jedi_traj_mod, only: fv3jedi_traj
use fv3jedi_traj_interface_mod, only: fv3jedi_traj_registry
use fv3jedi_geom_mod, only: fv3jedi_geom, fv3jedi_geom_registry
use fv3jedi_increment_mod, only: fv3jedi_increment, fv3jedi_increment_registry

implicit none
private

public :: fv3jedi_tlm_registry

! ------------------------------------------------------------------------------

!> Linked list interface - defines registry_t type
#define LISTED_TYPE fv3jedi_tlm
#include "linkedList_i.f"
type(registry_t) :: fv3jedi_tlm_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_create(c_conf, c_key_geom, c_key_self) bind (c,name='fv3jedi_tlm_create_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Key to tlm data
integer(c_int), intent(in)    :: c_key_geom  !< Geometry
type(c_ptr),    intent(in)    :: c_conf      !< pointer to object of class Config

type(fv3jedi_tlm),  pointer :: self
type(fv3jedi_geom), pointer :: geom

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_tlm_registry%init()
call fv3jedi_tlm_registry%add(c_key_self)
call fv3jedi_tlm_registry%get(c_key_self, self)

call tlm_create(self, geom, c_conf)

end subroutine c_fv3jedi_tlm_create

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_delete(c_key_self) bind (c,name='fv3jedi_tlm_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(fv3jedi_tlm), pointer :: self

call fv3jedi_tlm_registry%get(c_key_self, self)

call tlm_delete(self)

call fv3jedi_tlm_registry%remove(c_key_self)

end subroutine c_fv3jedi_tlm_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_initialize_tl(c_key_self, c_key_incr) bind(c,name='fv3jedi_tlm_initialize_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_self !< TLM
integer(c_int), intent(in) :: c_key_incr !< Increment

type(fv3jedi_tlm),       pointer :: self
type(fv3jedi_increment), pointer :: incr

call fv3jedi_increment_registry%get(c_key_incr,incr)
call fv3jedi_tlm_registry%get(c_key_self, self)

call tlm_initialize_tl(self, incr)

end subroutine c_fv3jedi_tlm_initialize_tl

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_step_tl(c_key_geom, c_key_self, c_key_incr, c_key_traj) bind(c,name='fv3jedi_tlm_step_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_self !< TLM
integer(c_int), intent(in) :: c_key_incr !< Increment
integer(c_int), intent(in) :: c_key_geom !< Geometry
integer(c_int), intent(in) :: c_key_traj !< Trajectory

type(fv3jedi_geom),      pointer :: geom
type(fv3jedi_tlm),       pointer :: self
type(fv3jedi_increment), pointer :: incr
type(fv3jedi_traj),      pointer :: traj

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_tlm_registry%get(c_key_self, self)
call fv3jedi_increment_registry%get(c_key_incr,incr)
call fv3jedi_traj_registry%get(c_key_traj,traj)

call tlm_step_tl(geom, self, incr, traj)

end subroutine c_fv3jedi_tlm_step_tl

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_finalize_tl(c_key_self, c_key_incr) bind(c,name='fv3jedi_tlm_finalize_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_self !< TLM
integer(c_int), intent(in) :: c_key_incr !< Increment

type(fv3jedi_tlm),       pointer :: self
type(fv3jedi_increment), pointer :: incr

call fv3jedi_increment_registry%get(c_key_incr,incr)
call fv3jedi_tlm_registry%get(c_key_self, self)

call tlm_finalize_tl(self, incr)

end subroutine c_fv3jedi_tlm_finalize_tl

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_initialize_ad(c_key_self, c_key_incr) bind(c,name='fv3jedi_tlm_initialize_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self !< TLM
integer(c_int), intent(in) :: c_key_incr !< Increment

type(fv3jedi_tlm),       pointer :: self
type(fv3jedi_increment), pointer :: incr

call fv3jedi_increment_registry%get(c_key_incr,incr)
call fv3jedi_tlm_registry%get(c_key_self, self)

call tlm_initialize_ad(self, incr)

end subroutine c_fv3jedi_tlm_initialize_ad

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_step_ad(c_key_geom, c_key_self, c_key_incr, c_key_traj) bind(c,name='fv3jedi_tlm_step_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self !< TLM
integer(c_int), intent(in) :: c_key_incr !< Increment
integer(c_int), intent(in) :: c_key_geom !< Geometry
integer(c_int), intent(in) :: c_key_traj !< Trajectory

type(fv3jedi_geom),      pointer :: geom
type(fv3jedi_tlm),       pointer :: self
type(fv3jedi_increment), pointer :: incr
type(fv3jedi_traj),      pointer :: traj

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_tlm_registry%get(c_key_self, self)
call fv3jedi_increment_registry%get(c_key_incr,incr)
call fv3jedi_traj_registry%get(c_key_traj,traj)

call tlm_step_ad(geom, self, incr, traj)

end subroutine c_fv3jedi_tlm_step_ad

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_finalize_ad(c_key_self, c_key_incr) bind(c,name='fv3jedi_tlm_finalize_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self !< TLM
integer(c_int), intent(in) :: c_key_incr !< Increment

type(fv3jedi_tlm),       pointer :: self
type(fv3jedi_increment), pointer :: incr

call fv3jedi_increment_registry%get(c_key_incr,incr)
call fv3jedi_tlm_registry%get(c_key_self, self)

call tlm_finalize_ad(self, incr)

end subroutine c_fv3jedi_tlm_finalize_ad

! ------------------------------------------------------------------------------

end module fv3jedi_tlm_interface_mod
