! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_tlm_interface_mod

! fckit
use fckit_configuration_module, only: fckit_configuration

! oops
use duration_mod
use iso_c_binding
use oops_variables_mod

! fv3-jedi
use fv3jedi_tlm_mod, only: fv3jedi_tlm
use fv3jedi_traj_mod, only: fv3jedi_traj
use fv3jedi_traj_interface_mod, only: fv3jedi_traj_registry
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_geom_interface_mod, only: fv3jedi_geom_registry
use fv3jedi_increment_mod, only: fv3jedi_increment
use fv3jedi_increment_interface_mod, only: fv3jedi_increment_registry

implicit none
private
public :: fv3jedi_tlm_registry

! --------------------------------------------------------------------------------------------------

!> Linked list interface - defines registry_t type
#define LISTED_TYPE fv3jedi_tlm
#include "oops/util/linkedList_i.f"
type(registry_t) :: fv3jedi_tlm_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_create(c_key_self, c_key_geom, c_conf) &
           bind (c,name='fv3jedi_tlm_create_f90')

implicit none
integer(c_int), intent(inout)  :: c_key_self  !< Key to tlm data
integer(c_int), intent(in)     :: c_key_geom  !< Geometry
type(c_ptr), value, intent(in) :: c_conf      !< pointer to object of class Config

type(fv3jedi_tlm),  pointer :: self
type(fv3jedi_geom), pointer :: geom
type(fckit_configuration)   :: f_conf

! Linked list
call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_tlm_registry%init()
call fv3jedi_tlm_registry%add(c_key_self)
call fv3jedi_tlm_registry%get(c_key_self, self)

! Fortran configuration
f_conf = fckit_configuration(c_conf)

! Implementation
call self%create(geom, f_conf)

end subroutine c_fv3jedi_tlm_create

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_delete(c_key_self) bind (c,name='fv3jedi_tlm_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self

type(fv3jedi_tlm), pointer :: self

! Linked list
call fv3jedi_tlm_registry%get(c_key_self, self)

! Implementation
call self%delete()

! Linked list
call fv3jedi_tlm_registry%remove(c_key_self)

end subroutine c_fv3jedi_tlm_delete

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_initialize_tl(c_key_self, c_key_incr, c_key_traj) &
           bind(c,name='fv3jedi_tlm_initialize_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_self !< TLM
integer(c_int), intent(in) :: c_key_incr !< Increment
integer(c_int), intent(in) :: c_key_traj !< Trajectory

type(fv3jedi_tlm),       pointer :: self
type(fv3jedi_increment), pointer :: incr
type(fv3jedi_traj),      pointer :: traj

! Linked list
call fv3jedi_increment_registry%get(c_key_incr,incr)
call fv3jedi_tlm_registry%get(c_key_self, self)
call fv3jedi_traj_registry%get(c_key_traj,traj)

! Implementation
call self%initialize_tl(incr, traj)

end subroutine c_fv3jedi_tlm_initialize_tl

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_step_tl(c_key_self, c_key_incr, c_key_traj) &
           bind(c,name='fv3jedi_tlm_step_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_self !< TLM
integer(c_int), intent(in) :: c_key_incr !< Increment
integer(c_int), intent(in) :: c_key_traj !< Trajectory

type(fv3jedi_tlm),       pointer :: self
type(fv3jedi_increment), pointer :: incr
type(fv3jedi_traj),      pointer :: traj

! Linked list
call fv3jedi_tlm_registry%get(c_key_self, self)
call fv3jedi_increment_registry%get(c_key_incr,incr)
call fv3jedi_traj_registry%get(c_key_traj,traj)

! Implementation
call self%step_tl(incr, traj)

end subroutine c_fv3jedi_tlm_step_tl

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_finalize_tl(c_key_self, c_key_incr) &
           bind(c,name='fv3jedi_tlm_finalize_tl_f90')

implicit none
integer(c_int), intent(in) :: c_key_self !< TLM
integer(c_int), intent(in) :: c_key_incr !< Increment

type(fv3jedi_tlm),       pointer :: self
type(fv3jedi_increment), pointer :: incr

! Linked list
call fv3jedi_increment_registry%get(c_key_incr,incr)
call fv3jedi_tlm_registry%get(c_key_self, self)

! Implementation
call self%finalize_tl(incr)

end subroutine c_fv3jedi_tlm_finalize_tl

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_initialize_ad(c_key_self, c_key_incr, c_key_traj) &
           bind(c,name='fv3jedi_tlm_initialize_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self !< TLM
integer(c_int), intent(in) :: c_key_incr !< Increment
integer(c_int), intent(in) :: c_key_traj !< Trajectory

type(fv3jedi_tlm),       pointer :: self
type(fv3jedi_increment), pointer :: incr
type(fv3jedi_traj),      pointer :: traj

! Linked list
call fv3jedi_increment_registry%get(c_key_incr,incr)
call fv3jedi_tlm_registry%get(c_key_self, self)
call fv3jedi_traj_registry%get(c_key_traj,traj)

! Implementation
call self%initialize_ad(incr, traj)

end subroutine c_fv3jedi_tlm_initialize_ad

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_step_ad(c_key_self, c_key_incr, c_key_traj) &
           bind(c,name='fv3jedi_tlm_step_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self !< TLM
integer(c_int), intent(in) :: c_key_incr !< Increment
integer(c_int), intent(in) :: c_key_traj !< Trajectory

type(fv3jedi_tlm),       pointer :: self
type(fv3jedi_increment), pointer :: incr
type(fv3jedi_traj),      pointer :: traj

! Linked list
call fv3jedi_tlm_registry%get(c_key_self, self)
call fv3jedi_increment_registry%get(c_key_incr,incr)
call fv3jedi_traj_registry%get(c_key_traj,traj)

! Implementation
call self%step_ad(incr, traj)

end subroutine c_fv3jedi_tlm_step_ad

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_tlm_finalize_ad(c_key_self, c_key_incr) &
           bind(c,name='fv3jedi_tlm_finalize_ad_f90')

implicit none
integer(c_int), intent(in) :: c_key_self !< TLM
integer(c_int), intent(in) :: c_key_incr !< Increment

type(fv3jedi_tlm),       pointer :: self
type(fv3jedi_increment), pointer :: incr

! Linked list
call fv3jedi_increment_registry%get(c_key_incr,incr)
call fv3jedi_tlm_registry%get(c_key_self, self)

! Implementation
call self%finalize_ad(incr)

end subroutine c_fv3jedi_tlm_finalize_ad

! --------------------------------------------------------------------------------------------------

end module fv3jedi_tlm_interface_mod
