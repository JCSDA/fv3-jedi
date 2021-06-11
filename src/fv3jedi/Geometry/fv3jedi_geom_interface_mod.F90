! (C) Copyright 2017-2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module handling geometry for the FV3 model

module fv3jedi_geom_interface_mod

use atlas_module, only: atlas_fieldset, atlas_functionspace
use iso_c_binding

use fckit_mpi_module,           only: fckit_mpi_comm
use fckit_configuration_module, only: fckit_configuration

use fields_metadata_mod, only: fields_metadata

use fv3jedi_kinds_mod
use fv3jedi_geom_mod

implicit none

private
public :: fv3jedi_geom_registry

! --------------------------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_geom

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_geom_registry

! --------------------------------------------------------------------------------------------------

contains

! --------------------------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_initialize(c_conf, c_comm) bind(c,name='fv3jedi_geom_initialize_f90')

implicit none

type(c_ptr), intent(in)        :: c_conf
type(c_ptr), value, intent(in) :: c_comm

type(fckit_mpi_comm)        :: f_comm
type(fckit_configuration)   :: f_conf

! Fortran APIs
! ------------
f_conf = fckit_configuration(c_conf)
f_comm = fckit_mpi_comm(c_comm)

call initialize(f_conf, f_comm)

end subroutine c_fv3jedi_geom_initialize

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_setup(c_key_self, c_conf, c_comm, c_fields_meta) &
                               bind(c, name='fv3jedi_geom_setup_f90')

implicit none

!Arguments
integer(c_int),     intent(inout) :: c_key_self
type(c_ptr),        intent(in)    :: c_conf
type(c_ptr), value, intent(in)    :: c_comm
type(c_ptr), value, intent(in)    :: c_fields_meta

type(fv3jedi_geom), pointer :: self
type(fckit_configuration)   :: f_conf
type(fckit_mpi_comm)        :: f_comm
type(fields_metadata)       :: f_fields_metadata

! LinkedList
! ----------
call fv3jedi_geom_registry%init()
call fv3jedi_geom_registry%add(c_key_self)
call fv3jedi_geom_registry%get(c_key_self,self)

! Fortran APIs
! ------------
f_conf            = fckit_configuration(c_conf)
f_comm            = fckit_mpi_comm(c_comm)
f_fields_metadata = fields_metadata(c_fields_meta)

! Call implementation
! -------------------
call self%create(f_conf, f_comm, f_fields_metadata)

end subroutine c_fv3jedi_geom_setup

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_clone(c_key_self, c_key_other, c_fields_meta) bind(c,name='fv3jedi_geom_clone_f90')

implicit none

integer(c_int),     intent(inout) :: c_key_self
integer(c_int),     intent(in)    :: c_key_other
type(c_ptr), value, intent(in)    :: c_fields_meta

type(fv3jedi_geom), pointer :: self
type(fv3jedi_geom), pointer :: other
type(fields_metadata)       :: f_fields_metadata

! LinkedList
! ----------
call fv3jedi_geom_registry%add(c_key_self)
call fv3jedi_geom_registry%get(c_key_self, self)
call fv3jedi_geom_registry%get(c_key_other, other)

! Fortran APIs
! ------------
f_fields_metadata = fields_metadata(c_fields_meta)

! Call implementation
! -------------------
call self%clone(other, f_fields_metadata)

end subroutine c_fv3jedi_geom_clone

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_delete(c_key_self) bind(c,name='fv3jedi_geom_delete_f90')

implicit none

integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_geom), pointer :: self

! LinkedList
! ----------
call fv3jedi_geom_registry%get(c_key_self, self)

! Call implementation
! -------------------
call self%delete()

! LinkedList
! ----------
call fv3jedi_geom_registry%remove(c_key_self)

end subroutine c_fv3jedi_geom_delete

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_print(c_key_self, c_cube) bind(c,name='fv3jedi_geom_print_f90')

implicit none

integer(c_int), intent(in)    :: c_key_self
integer(c_int), intent(inout) :: c_cube

type(fv3jedi_geom), pointer :: self

! LinkedList
! ----------
call fv3jedi_geom_registry%get(c_key_self, self)

! Get Cube size
! -------------
c_cube = self%npx - 1

end subroutine c_fv3jedi_geom_print

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_set_atlas_lonlat(c_key_self, c_afieldset) &
                                           bind(c,name='fv3jedi_geom_set_atlas_lonlat_f90')

implicit none

!Arguments
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in), value :: c_afieldset

type(fv3jedi_geom), pointer :: self
type(atlas_fieldset) :: afieldset

! LinkedList
! ----------
call fv3jedi_geom_registry%get(c_key_self,self)

! Fortran APIs
! ------------
afieldset = atlas_fieldset(c_afieldset)

! Call implementation
! -------------------
call self%set_atlas_lonlat(afieldset)

end subroutine c_fv3jedi_geom_set_atlas_lonlat

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_set_atlas_functionspace_pointer(c_key_self,c_afunctionspace) &
                                                          bind(c,name='fv3jedi_geom_set_atlas_functionspace_pointer_f90')

integer(c_int), intent(in)     :: c_key_self
type(c_ptr), intent(in), value :: c_afunctionspace

type(fv3jedi_geom),pointer :: self

! LinkedList
! ----------
call fv3jedi_geom_registry%get(c_key_self,self)

! Create function space
! ---------------------
self%afunctionspace = atlas_functionspace(c_afunctionspace)

end subroutine c_fv3jedi_geom_set_atlas_functionspace_pointer

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_fill_atlas_fieldset(c_key_self, c_afieldset) &
                                              bind(c,name='fv3jedi_geom_fill_atlas_fieldset_f90')

implicit none

integer(c_int),     intent(in) :: c_key_self
type(c_ptr), value, intent(in) :: c_afieldset

type(fv3jedi_geom), pointer :: self
type(atlas_fieldset) :: afieldset

! LinkedList
! ----------
call fv3jedi_geom_registry%get(c_key_self,self)
afieldset = atlas_fieldset(c_afieldset)

! Call implementation
! -------------------
call self%fill_atlas_fieldset(afieldset)

end subroutine c_fv3jedi_geom_fill_atlas_fieldset

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_start_end(c_key_self, ist, iend, jst, jend, npz) &
                                    bind(c, name='fv3jedi_geom_start_end_f90')

implicit none

integer(c_int), intent( in) :: c_key_self
integer(c_int), intent(out) :: ist, iend, jst, jend, npz

type(fv3jedi_geom), pointer :: self

! LinkedList
! ----------
call fv3jedi_geom_registry%get(c_key_self, self)

ist  = self%isc
iend = self%iec
jst  = self%jsc
jend = self%jec
npz  = self%npz

end subroutine c_fv3jedi_geom_start_end

!--------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_verticalCoord(c_key_self, vc, npz, psurf) &
                                    bind(c, name='fv3jedi_geom_verticalCoord_f90')

implicit none

integer(c_int),    intent( in) :: c_key_self
integer(c_int),    intent( in) :: npz
real(c_double), intent( in) :: psurf
real(c_double), intent(out) :: vc(npz)

type(fv3jedi_geom), pointer :: self

! LinkedList
! ----------
call fv3jedi_geom_registry%get(c_key_self, self)

! Call implementation
! -------------------
call getVerticalCoordLogP(self, vc, npz, psurf)

end subroutine c_fv3jedi_geom_verticalCoord

! --------------------------------------------------------------------------------------------------

end module fv3jedi_geom_interface_mod
