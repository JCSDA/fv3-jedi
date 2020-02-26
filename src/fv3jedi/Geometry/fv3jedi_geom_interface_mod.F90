! (C) Copyright 2017-2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Fortran module handling geometry for the FV3 model

module fv3jedi_geom_interface_mod

use atlas_module
use fv3jedi_kinds_mod
use iso_c_binding
use fv3jedi_geom_mod

implicit none
private

public :: fv3jedi_geom_registry

! ------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_geom

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_geom_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geo_setup(c_key_self, c_conf, c_comm) bind(c,name='fv3jedi_geo_setup_f90')
use fckit_mpi_module,   only: fckit_mpi_comm

implicit none

!Arguments
integer(c_int), intent(inout)  :: c_key_self
type(c_ptr), intent(in)        :: c_conf
type(c_ptr), value, intent(in) :: c_comm

type(fv3jedi_geom), pointer :: self
type(fckit_mpi_comm)        :: f_comm

! Init, add and get key
! ---------------------
call fv3jedi_geom_registry%init()
call fv3jedi_geom_registry%add(c_key_self)
call fv3jedi_geom_registry%get(c_key_self,self)

f_comm = fckit_mpi_comm(c_comm)
call create(self,c_conf, f_comm)

end subroutine c_fv3jedi_geo_setup

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geo_create_atlas_grid_conf(c_key_self, c_conf) bind(c,name='fv3jedi_geo_create_atlas_grid_conf_f90')

implicit none

!Arguments
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in)    :: c_conf

type(fv3jedi_geom), pointer :: self

! Get key
call fv3jedi_geom_registry%get(c_key_self,self)

call create_atlas_grid_conf(self,c_conf)

end subroutine c_fv3jedi_geo_create_atlas_grid_conf

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geo_set_atlas_functionspace_pointer(c_key_self,c_afunctionspace) &
 & bind(c,name='fv3jedi_geo_set_atlas_functionspace_pointer_f90')

!Arguments
integer(c_int), intent(in)     :: c_key_self
type(c_ptr), intent(in), value :: c_afunctionspace

type(fv3jedi_geom),pointer :: self

! Get key
call fv3jedi_geom_registry%get(c_key_self,self)
self%afunctionspace = atlas_functionspace_nodecolumns(c_afunctionspace)

end subroutine c_fv3jedi_geo_set_atlas_functionspace_pointer

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geo_fill_atlas_fieldset(c_key_self, c_afieldset) &
 & bind(c,name='fv3jedi_geo_fill_atlas_fieldset_f90')

implicit none

!Arguments
integer(c_int), intent(in)     :: c_key_self
type(c_ptr), intent(in), value :: c_afieldset

type(fv3jedi_geom), pointer :: self
type(atlas_fieldset) :: afieldset

! Get key
call fv3jedi_geom_registry%get(c_key_self,self)
afieldset = atlas_fieldset(c_afieldset)

call fill_atlas_fieldset(self,afieldset)

end subroutine c_fv3jedi_geo_fill_atlas_fieldset

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geo_clone(c_key_self, c_key_other) bind(c,name='fv3jedi_geo_clone_f90')

implicit none

integer(c_int), intent(in   ) :: c_key_self
integer(c_int), intent(inout) :: c_key_other

type(fv3jedi_geom), pointer :: self, other

!add, get, get key
call fv3jedi_geom_registry%add(c_key_other)
call fv3jedi_geom_registry%get(c_key_other, other)
call fv3jedi_geom_registry%get(c_key_self, self)

call clone(self, other)

end subroutine c_fv3jedi_geo_clone

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geo_delete(c_key_self) bind(c,name='fv3jedi_geo_delete_f90')

implicit none

integer(c_int), intent(inout) :: c_key_self
type(fv3jedi_geom), pointer :: self

! Get key
call fv3jedi_geom_registry%get(c_key_self, self)

call delete(self)

! Remove key
call fv3jedi_geom_registry%remove(c_key_self)

end subroutine c_fv3jedi_geo_delete

! ------------------------------------------------------------------------------

subroutine c_fv3jedi_geo_info(c_key_self) bind(c,name='fv3jedi_geo_info_f90')

implicit none

integer(c_int), intent(in   ) :: c_key_self
type(fv3jedi_geom), pointer :: self

call fv3jedi_geom_registry%get(c_key_self, self)

call info(self)

end subroutine c_fv3jedi_geo_info

! ------------------------------------------------------------------------------
!> return begin and end of local geometry
subroutine c_fv3jedi_geo_start_end(c_key_self, ist, iend, jst, jend, npz) bind(c, name='fv3jedi_geo_start_end_f90')

  implicit none

  integer(c_int), intent( in) :: c_key_self
  integer(c_int), intent(out) :: ist, iend, jst, jend, npz

  type(fv3jedi_geom), pointer :: self
  call fv3jedi_geom_registry%get(c_key_self, self)

  ist  = self%isc
  iend = self%iec
  jst  = self%jsc
  jend = self%jec
  npz  = self%npz

end subroutine c_fv3jedi_geo_start_end



! ------------------------------------------------------------------------------

end module fv3jedi_geom_interface_mod
