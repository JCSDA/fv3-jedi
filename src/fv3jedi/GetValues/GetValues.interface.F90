! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! --------------------------------------------------------------------------------------------------

module fv3jedi_getvalues_interface_mod

! Intrinsic
use iso_c_binding

! oops dependencies
use datetime_mod
use duration_mod
use oops_variables_mod

! ufo dependencies
use ufo_locs_mod
use ufo_locs_mod_c, only: ufo_locs_registry
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry

! self dependency
use fv3jedi_getvalues_mod, only: fv3jedi_getvalues

! fv3jedi dependencies
use fv3jedi_geom_interface_mod,  only: fv3jedi_geom_registry
use fv3jedi_geom_mod,            only: fv3jedi_geom
use fv3jedi_kinds_mod,           only: kind_real
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry
use fv3jedi_state_mod,           only: fv3jedi_state

implicit none
private
public :: fv3jedi_getvalues_registry

! --------------------------------------------------------------------------------------------------

!> Linked list interface
#define LISTED_TYPE fv3jedi_getvalues
#include "oops/util/linkedList_i.f"
type(registry_t) :: fv3jedi_getvalues_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_getvalues_create_c(c_key_self, c_key_geom, c_key_locs) &
           bind (c, name='fv3jedi_getvalues_create_f90')

integer(c_int),     intent(inout) :: c_key_self      !< Key to self
integer(c_int),     intent(in)    :: c_key_geom      !< Key to geometry
integer(c_int),     intent(in)    :: c_key_locs      !< Key to observation locations

type(fv3jedi_getvalues), pointer :: self
type(fv3jedi_geom),      pointer :: geom
type(ufo_locs),          pointer :: locs

! Create object
call fv3jedi_getvalues_registry%init()
call fv3jedi_getvalues_registry%add(c_key_self)
call fv3jedi_getvalues_registry%get(c_key_self, self)

! Others
call fv3jedi_geom_registry%get(c_key_geom, geom)
call ufo_locs_registry%get(c_key_locs, locs)

! Call method
call self%create(geom, locs)

end subroutine fv3jedi_getvalues_create_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_getvalues_delete_c(c_key_self) bind (c, name='fv3jedi_getvalues_delete_f90')

integer(c_int), intent(inout) :: c_key_self !< Key to self

type(fv3jedi_getvalues), pointer :: self

! Get object
call fv3jedi_getvalues_registry%get(c_key_self, self)

! Call method
call self%delete()

! Remove object
call fv3jedi_getvalues_registry%remove(c_key_self)

end subroutine fv3jedi_getvalues_delete_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_getvalues_fill_geovals_c(c_key_self, c_key_geom, c_key_state, c_t1, c_t2, &
                                            c_key_locs, c_key_geovals) &
           bind (c, name='fv3jedi_getvalues_fill_geovals_f90')

integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_state
type(c_ptr),    intent(in) :: c_t1
type(c_ptr),    intent(in) :: c_t2
integer(c_int), intent(in) :: c_key_locs
integer(c_int), intent(in) :: c_key_geovals

type(fv3jedi_getvalues), pointer :: self
type(fv3jedi_geom),      pointer :: geom
type(fv3jedi_state),     pointer :: state
type(datetime)                   :: t1
type(datetime)                   :: t2
type(ufo_locs),          pointer :: locs
type(ufo_geovals),       pointer :: geovals

! Get objects
call fv3jedi_getvalues_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_state_registry%get(c_key_state, state)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
call ufo_locs_registry%get(c_key_locs, locs)
call ufo_geovals_registry%get(c_key_geovals, geovals)

! Call method
call self%fill_geovals(geom, state%fields, t1, t2, locs, geovals)

end subroutine fv3jedi_getvalues_fill_geovals_c

! --------------------------------------------------------------------------------------------------

end module fv3jedi_getvalues_interface_mod
