! (C) Copyright 2021-2022 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! --------------------------------------------------------------------------------------------------

module fv3jedi_io_cube_sphere_history_interface_mod

! iso
use iso_c_binding

! fckit
use fckit_configuration_module,      only: fckit_configuration

! fv3-jedi
use fv3jedi_increment_mod,              only: fv3jedi_increment
use fv3jedi_increment_interface_mod,    only: fv3jedi_increment_registry
use fv3jedi_geom_mod,                   only: fv3jedi_geom
use fv3jedi_geom_interface_mod,         only: fv3jedi_geom_registry
use fv3jedi_io_cube_sphere_history_mod, only: fv3jedi_io_cube_sphere_history
use fv3jedi_state_mod,                  only: fv3jedi_state
use fv3jedi_state_interface_mod,        only: fv3jedi_state_registry
use fv3jedi_geom_mod,                   only: fv3jedi_geom
use fv3jedi_geom_interface_mod,         only: fv3jedi_geom_registry

implicit none
private
public :: fv3jedi_io_cube_sphere_history_registry

! --------------------------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_io_cube_sphere_history

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_io_cube_sphere_history_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_io_cube_sphere_history_create(c_key_self, c_conf, c_key_geom) &
           bind (c,name='fv3jedi_io_cube_sphere_history_create_f90')

integer(c_int),     intent(inout) :: c_key_self
type(c_ptr), value, intent(in)    :: c_conf
integer(c_int),     intent(in)    :: c_key_geom

type(fv3jedi_io_cube_sphere_history), pointer :: f_self
type(fckit_configuration)                     :: f_conf
type(fv3jedi_geom),                   pointer :: f_geom

! Linked list
! -----------
call fv3jedi_io_cube_sphere_history_registry%init()
call fv3jedi_io_cube_sphere_history_registry%add(c_key_self)
call fv3jedi_io_cube_sphere_history_registry%get(c_key_self, f_self)

call fv3jedi_geom_registry%get(c_key_geom, f_geom)

! Fortran APIs
! ------------
f_conf = fckit_configuration(c_conf)

! Call implementation
! -------------------
call f_self%create(f_geom, f_conf)

end subroutine c_fv3jedi_io_cube_sphere_history_create

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_io_cube_sphere_history_delete(c_key_self) &
           bind (c,name='fv3jedi_io_cube_sphere_history_delete_f90')

integer(c_int), intent(inout) :: c_key_self  !< Change variable structure

type(fv3jedi_io_cube_sphere_history), pointer :: f_self

! Linked list
! -----------
call fv3jedi_io_cube_sphere_history_registry%get(c_key_self, f_self)

! Call implementation
! -------------------
call f_self%delete()

! Linked list
! -----------
call fv3jedi_io_cube_sphere_history_registry%remove(c_key_self)

end subroutine c_fv3jedi_io_cube_sphere_history_delete

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_io_cube_sphere_history_read_state(c_key_self, c_key_state) &
           bind (c,name='fv3jedi_io_cube_sphere_history_read_state_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_state

type(fv3jedi_io_cube_sphere_history), pointer :: f_self
type(fv3jedi_state),                  pointer :: f_state

! Linked list
! -----------
call fv3jedi_io_cube_sphere_history_registry%get(c_key_self, f_self)
call fv3jedi_state_registry%get(c_key_state, f_state)

! Call implementation
! -------------------
call f_self%read(f_state%time, f_state%fields)

end subroutine c_fv3jedi_io_cube_sphere_history_read_state

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_io_cube_sphere_history_read_increment(c_key_self, c_key_increment) &
           bind (c,name='fv3jedi_io_cube_sphere_history_read_increment_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_increment

type(fv3jedi_io_cube_sphere_history),   pointer :: f_self
type(fv3jedi_increment),                pointer :: f_increment

! Linked list
! -----------
call fv3jedi_io_cube_sphere_history_registry%get(c_key_self, f_self)
call fv3jedi_increment_registry%get(c_key_increment, f_increment)

! Call implementation
! -------------------
call f_self%read(f_increment%time, f_increment%fields)

end subroutine c_fv3jedi_io_cube_sphere_history_read_increment

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_io_cube_sphere_history_write_state(c_key_self, c_key_state) &
           bind (c,name='fv3jedi_io_cube_sphere_history_write_state_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_state

type(fv3jedi_io_cube_sphere_history), pointer :: f_self
type(fv3jedi_state),                  pointer :: f_state

! Linked list
! -----------
call fv3jedi_io_cube_sphere_history_registry%get(c_key_self, f_self)
call fv3jedi_state_registry%get(c_key_state, f_state)

! Call implementation
! -------------------
call f_self%write(f_state%fields, f_state%time)

end subroutine c_fv3jedi_io_cube_sphere_history_write_state

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_io_cube_sphere_history_write_increment(c_key_self, c_key_increment) &
           bind (c,name='fv3jedi_io_cube_sphere_history_write_increment_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_increment

type(fv3jedi_io_cube_sphere_history),   pointer :: f_self
type(fv3jedi_increment),                pointer :: f_increment

! Linked list
! -----------
call fv3jedi_io_cube_sphere_history_registry%get(c_key_self, f_self)
call fv3jedi_increment_registry%get(c_key_increment, f_increment)

! Call implementation
! -------------------
call f_self%write(f_increment%fields, f_increment%time)

end subroutine c_fv3jedi_io_cube_sphere_history_write_increment

! --------------------------------------------------------------------------------------------------

end module fv3jedi_io_cube_sphere_history_interface_mod
