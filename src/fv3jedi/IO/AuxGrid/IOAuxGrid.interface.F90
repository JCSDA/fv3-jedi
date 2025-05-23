! (C) Copyright 2022 NOAA
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

! --------------------------------------------------------------------------------------------------

module fv3jedi_io_auxgrid_interface_mod

! iso
use iso_c_binding

! fckit
use fckit_configuration_module,      only: fckit_configuration

! fv3-jedi
use fv3jedi_increment_mod,           only: fv3jedi_increment
use fv3jedi_increment_interface_mod, only: fv3jedi_increment_registry
use fv3jedi_geom_mod,                only: fv3jedi_geom
use fv3jedi_geom_interface_mod,      only: fv3jedi_geom_registry
use fv3jedi_io_auxgrid_mod,           only: fv3jedi_io_auxgrid
use fv3jedi_state_mod,               only: fv3jedi_state
use fv3jedi_state_interface_mod,     only: fv3jedi_state_registry
use fv3jedi_geom_mod,                only: fv3jedi_geom
use fv3jedi_geom_interface_mod,      only: fv3jedi_geom_registry

implicit none
private
public :: fv3jedi_io_auxgrid_registry

! --------------------------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_io_auxgrid

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_io_auxgrid_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_io_auxgrid_create(c_key_self, c_conf, c_key_geom) &
           bind (c,name='fv3jedi_io_auxgrid_create_f90')

integer(c_int),     intent(inout) :: c_key_self
type(c_ptr), value, intent(in)    :: c_conf
integer(c_int),     intent(in)    :: c_key_geom

type(fv3jedi_io_auxgrid), pointer :: f_self
type(fckit_configuration)        :: f_conf
type(fv3jedi_geom),      pointer :: f_geom

! Linked list
! -----------
call fv3jedi_io_auxgrid_registry%init()
call fv3jedi_io_auxgrid_registry%add(c_key_self)
call fv3jedi_io_auxgrid_registry%get(c_key_self, f_self)

call fv3jedi_geom_registry%get(c_key_geom, f_geom)

! Fortran APIs
! ------------
f_conf = fckit_configuration(c_conf)

! Call implementation
! -------------------
call f_self%create(f_geom, f_conf)

end subroutine c_fv3jedi_io_auxgrid_create

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_io_auxgrid_delete(c_key_self) bind (c,name='fv3jedi_io_auxgrid_delete_f90')

integer(c_int), intent(inout) :: c_key_self  !< Change variable structure

type(fv3jedi_io_auxgrid), pointer :: f_self

! Linked list
! -----------
call fv3jedi_io_auxgrid_registry%get(c_key_self, f_self)

! Call implementation
! -------------------
call f_self%delete()

! Linked list
! -----------
call fv3jedi_io_auxgrid_registry%remove(c_key_self)

end subroutine c_fv3jedi_io_auxgrid_delete

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_io_auxgrid_write_state(c_key_self, c_key_state, c_fileionamesconf, &
                                            c_fileioscalingconf) &
           bind (c,name='fv3jedi_io_auxgrid_write_state_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_state
type(c_ptr), value, intent(in) :: c_fileionamesconf
type(c_ptr), value, intent(in) :: c_fileioscalingconf

type(fv3jedi_io_auxgrid), pointer :: f_self
type(fv3jedi_state),      pointer :: f_state
type(fckit_configuration)         :: f_fileionamesconf
type(fckit_configuration)         :: f_fileioscalingconf

! Linked list
! -----------
call fv3jedi_io_auxgrid_registry%get(c_key_self, f_self)
call fv3jedi_state_registry%get(c_key_state, f_state)

! APIS
! ----
f_fileionamesconf = fckit_configuration(c_fileionamesconf)
f_fileioscalingconf = fckit_configuration(c_fileioscalingconf)

! Call implementation
! -------------------
call f_self%write(f_state%time, f_state%fields, f_fileionamesconf, f_fileioscalingconf)

end subroutine c_fv3jedi_io_auxgrid_write_state

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_io_auxgrid_write_increment(c_key_self, c_key_increment, c_fileionamesconf, &
                                                c_fileioscalingconf) &
           bind (c,name='fv3jedi_io_auxgrid_write_increment_f90')

implicit none
integer(c_int),     intent(in) :: c_key_self
integer(c_int),     intent(in) :: c_key_increment
type(c_ptr), value, intent(in) :: c_fileionamesconf
type(c_ptr), value, intent(in) :: c_fileioscalingconf

type(fv3jedi_io_auxgrid), pointer :: f_self
type(fv3jedi_increment), pointer  :: f_increment
type(fckit_configuration)         :: f_fileionamesconf
type(fckit_configuration)         :: f_fileioscalingconf

! Linked list
! -----------
call fv3jedi_io_auxgrid_registry%get(c_key_self, f_self)
call fv3jedi_increment_registry%get(c_key_increment, f_increment)

! APIS
! ----
f_fileionamesconf = fckit_configuration(c_fileionamesconf)
f_fileioscalingconf = fckit_configuration(c_fileioscalingconf)

! Call implementation
! -------------------
call f_self%write(f_increment%time, f_increment%fields, f_fileionamesconf, f_fileioscalingconf)

end subroutine c_fv3jedi_io_auxgrid_write_increment

! --------------------------------------------------------------------------------------------------

end module fv3jedi_io_auxgrid_interface_mod
