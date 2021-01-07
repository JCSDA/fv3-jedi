! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! --------------------------------------------------------------------------------------------------

module fv3jedi_lineargetvalues_interface_mod

! Intrinsic
use iso_c_binding

! oops dependencies
use datetime_mod

! ufo dependencies
use ufo_locations_mod
use ufo_geovals_mod
use ufo_geovals_mod_c, only: ufo_geovals_registry

! self dependency
use fv3jedi_lineargetvalues_mod, only: fv3jedi_lineargetvalues

! fv3jedi dependencies
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_geom_interface_mod, only: fv3jedi_geom_registry
use fv3jedi_increment_mod, only: fv3jedi_increment
use fv3jedi_increment_interface_mod, only: fv3jedi_increment_registry
use fv3jedi_state_mod, only: fv3jedi_state
use fv3jedi_state_interface_mod, only: fv3jedi_state_registry

implicit none
private
public :: fv3jedi_lineargetvalues_registry

! --------------------------------------------------------------------------------------------------

!> Linked list interface
#define LISTED_TYPE fv3jedi_lineargetvalues
#include "oops/util/linkedList_i.f"
type(registry_t) :: fv3jedi_lineargetvalues_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_lineargetvalues_create_c(c_key_self, c_key_geom, c_locs) &
           bind (c,name='fv3jedi_lineargetvalues_create_f90')
integer(c_int),     intent(inout) :: c_key_self      !< Key to self
integer(c_int),     intent(in)    :: c_key_geom      !< Key to geometry
type(c_ptr), value, intent(in)    :: c_locs          !< Observation locations

type(fv3jedi_lineargetvalues), pointer :: self
type(fv3jedi_geom),            pointer :: geom
type(ufo_locations)                    :: locs

! Create object
call fv3jedi_lineargetvalues_registry%init()
call fv3jedi_lineargetvalues_registry%add(c_key_self)
call fv3jedi_lineargetvalues_registry%get(c_key_self, self)

! Others
call fv3jedi_geom_registry%get(c_key_geom, geom)
locs = ufo_locations(c_locs)

! Call method
call self%create(geom, locs)

end subroutine fv3jedi_lineargetvalues_create_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_lineargetvalues_delete_c(c_key_self) &
           bind (c,name='fv3jedi_lineargetvalues_delete_f90')

integer(c_int), intent(inout) :: c_key_self !< Key to self

type(fv3jedi_lineargetvalues), pointer :: self

! Get object
call fv3jedi_lineargetvalues_registry%get(c_key_self, self)

! Call method
call self%delete()

! Remove object
call fv3jedi_lineargetvalues_registry%remove(c_key_self)

end subroutine fv3jedi_lineargetvalues_delete_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_lineargetvalues_set_trajectory_c(c_key_self, c_key_geom, c_key_state, c_t1, &
                                                    c_t2, c_locs, c_key_geovals) &
           bind (c,name='fv3jedi_lineargetvalues_set_trajectory_f90')

integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_state
type(c_ptr), value, intent(in) :: c_t1
type(c_ptr), value, intent(in) :: c_t2
type(c_ptr), value, intent(in) :: c_locs
integer(c_int), intent(in) :: c_key_geovals

type(fv3jedi_lineargetvalues), pointer :: self
type(fv3jedi_geom),            pointer :: geom
type(fv3jedi_state),           pointer :: state
type(datetime)                         :: t1
type(datetime)                         :: t2
type(ufo_locations)                    :: locs
type(ufo_geovals),             pointer :: geovals

! Get objects
call fv3jedi_lineargetvalues_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_state_registry%get(c_key_state, state)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
locs = ufo_locations(c_locs)
call ufo_geovals_registry%get(c_key_geovals, geovals)

! Call method
call self%set_trajectory(geom, state%fields, t1, t2, locs, geovals)

end subroutine fv3jedi_lineargetvalues_set_trajectory_c


! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_lineargetvalues_fill_geovals_tl_c(c_key_self, c_key_geom, c_key_inc, c_t1, &
                                                     c_t2, c_locs, c_key_geovals) &
           bind (c,name='fv3jedi_lineargetvalues_fill_geovals_tl_f90')

integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_inc
type(c_ptr), value, intent(in) :: c_t1
type(c_ptr), value, intent(in) :: c_t2
type(c_ptr), value, intent(in) :: c_locs
integer(c_int), intent(in) :: c_key_geovals

type(fv3jedi_lineargetvalues), pointer :: self
type(fv3jedi_geom),            pointer :: geom
type(fv3jedi_increment),       pointer :: inc
type(datetime)                         :: t1
type(datetime)                         :: t2
type(ufo_locations)                    :: locs
type(ufo_geovals),             pointer :: geovals

! Get objects
call fv3jedi_lineargetvalues_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_increment_registry%get(c_key_inc, inc)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
locs = ufo_locations(c_locs)
call ufo_geovals_registry%get(c_key_geovals, geovals)

! Call method
call self%fill_geovals_tl(geom, inc%fields, t1, t2, locs, geovals)

end subroutine fv3jedi_lineargetvalues_fill_geovals_tl_c

! --------------------------------------------------------------------------------------------------

subroutine fv3jedi_lineargetvalues_fill_geovals_ad_c(c_key_self, c_key_geom, c_key_inc, c_t1, &
                                                     c_t2, c_locs, c_key_geovals) &
           bind (c,name='fv3jedi_lineargetvalues_fill_geovals_ad_f90')

integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_geom
integer(c_int), intent(in) :: c_key_inc
type(c_ptr), value, intent(in) :: c_t1
type(c_ptr), value, intent(in) :: c_t2
type(c_ptr), value, intent(in) :: c_locs
integer(c_int), intent(in) :: c_key_geovals

type(fv3jedi_lineargetvalues), pointer :: self
type(fv3jedi_geom),            pointer :: geom
type(fv3jedi_increment),       pointer :: inc
type(datetime)                         :: t1
type(datetime)                         :: t2
type(ufo_locations)                    :: locs
type(ufo_geovals),             pointer :: geovals

! Get objects
call fv3jedi_lineargetvalues_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_geom, geom)
call fv3jedi_increment_registry%get(c_key_inc, inc)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
locs = ufo_locations(c_locs)
call ufo_geovals_registry%get(c_key_geovals, geovals)

! Call method
call self%fill_geovals_ad(geom, inc%fields, t1, t2, locs, geovals)

end subroutine fv3jedi_lineargetvalues_fill_geovals_ad_c

! --------------------------------------------------------------------------------------------------

end module fv3jedi_lineargetvalues_interface_mod
