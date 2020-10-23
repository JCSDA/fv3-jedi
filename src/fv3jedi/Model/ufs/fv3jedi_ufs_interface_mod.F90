! (C) Copyright 2020 NOAA
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module fv3jedi_ufs_interface_mod

! iso
use iso_c_binding

! oops
use datetime_mod
use duration_mod

! fckit
use fckit_configuration_module,  only: fckit_configuration

! fv3jedi
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_geom_interface_mod,  only: fv3jedi_geom_registry
use fv3jedi_state_mod,           only: fv3jedi_state
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry
use fv3jedi_ufs_mod,             only: model_ufs

implicit none
private
public :: fv3jedi_ufs_registry

!> Linked list interface
#define LISTED_TYPE model_ufs
#include "oops/util/linkedList_i.f"
type(registry_t) :: fv3jedi_ufs_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_ufs_create(c_key_self, c_conf, c_key_geom) &
           bind (c,name='fv3jedi_ufs_create_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self  !< Key to model data
integer(c_int), intent(in)    :: c_key_geom  !< Geometry
type(c_ptr),    intent(in)    :: c_conf      !< pointer to object of class Config

type(model_ufs),     pointer :: self
type(fv3jedi_geom),  pointer :: geom
type(fckit_configuration)    :: f_conf

call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_ufs_registry%init()
call fv3jedi_ufs_registry%add(c_key_self)
call fv3jedi_ufs_registry%get(c_key_self, self)

f_conf = fckit_configuration(c_conf)

call self%create(f_conf, geom)

end subroutine c_fv3jedi_ufs_create

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_ufs_delete(c_key_self) bind (c,name='fv3jedi_ufs_delete_f90')

implicit none
integer(c_int), intent(inout) :: c_key_self
type(model_ufs), pointer :: self

call fv3jedi_ufs_registry%get(c_key_self, self)

call self%delete()

call fv3jedi_ufs_registry%remove(c_key_self)

end subroutine c_fv3jedi_ufs_delete

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_ufs_initialize(c_key_self, c_key_state, c_dt) &
           bind(c,name='fv3jedi_ufs_initialize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state
type(c_ptr),    intent(in) :: c_dt        !< DateTime

type(model_ufs),     pointer :: self
type(fv3jedi_state), pointer :: state
type(datetime) :: fdate

call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_ufs_registry%get(c_key_self, self)

call c_f_datetime(c_dt, fdate)

call self%initialize(state, fdate)

end subroutine c_fv3jedi_ufs_initialize

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_ufs_step(c_key_self, c_key_state, c_dt1, c_dt2) &
           bind(c,name='fv3jedi_ufs_step_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state
type(c_ptr),    intent(in) :: c_dt1       !< DateTime
type(c_ptr),    intent(in) :: c_dt2       !< DateTime

type(model_ufs), pointer :: self
type(fv3jedi_state), pointer :: state

type(datetime) :: fdate1
type(datetime) :: fdate2
type(duration) :: dt
integer(c_int) :: ts
character(len=20) :: vdatestrz

call fv3jedi_ufs_registry%get(c_key_self, self)
call fv3jedi_state_registry%get(c_key_state,state)

call c_f_datetime(c_dt1, fdate1)
call datetime_to_string(fdate1, vdatestrz)
call c_f_datetime(c_dt2, fdate2)
call datetime_to_string(fdate2, vdatestrz)
call self%step(state, fdate1, fdate2)

end subroutine c_fv3jedi_ufs_step

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_ufs_finalize(c_key_self, c_key_state, c_dt) &
           bind(c,name='fv3jedi_ufs_finalize_f90')

implicit none
integer(c_int), intent(in) :: c_key_self  !< Model
integer(c_int), intent(in) :: c_key_state !< Model state
type(c_ptr),    intent(in) :: c_dt        !< DateTime

type(model_ufs),     pointer :: self
type(fv3jedi_state), pointer :: state

type(datetime) :: fdate

call fv3jedi_state_registry%get(c_key_state,state)
call fv3jedi_ufs_registry%get(c_key_self, self)

call c_f_datetime(c_dt, fdate)

call self%finalize(state, fdate)

end subroutine c_fv3jedi_ufs_finalize

! --------------------------------------------------------------------------------------------------

end module fv3jedi_ufs_interface_mod
