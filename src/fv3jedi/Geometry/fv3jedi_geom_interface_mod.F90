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

type(c_ptr), value, intent(in) :: c_conf
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

subroutine c_fv3jedi_geom_setup(c_key_self, c_conf, c_comm, c_nlev) &
                               bind(c, name='fv3jedi_geom_setup_f90')

!Arguments
integer(c_int),     intent(inout) :: c_key_self
type(c_ptr), value, intent(in)    :: c_conf
type(c_ptr), value, intent(in)    :: c_comm
integer(c_int),     intent(inout) :: c_nlev

type(fv3jedi_geom), pointer :: self
type(fckit_configuration)   :: f_conf
type(fckit_mpi_comm)        :: f_comm
integer                     :: f_nlev

! LinkedList
! ----------
call fv3jedi_geom_registry%init()
call fv3jedi_geom_registry%add(c_key_self)
call fv3jedi_geom_registry%get(c_key_self,self)

! Fortran APIs
! ------------
f_conf            = fckit_configuration(c_conf)
f_comm            = fckit_mpi_comm(c_comm)

! Call implementation
! -------------------
call self%create(f_conf, f_comm, f_nlev)

! Pass number of levels
! ---------------------
c_nlev = f_nlev

end subroutine c_fv3jedi_geom_setup

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_addfmd_setup(c_key_self, c_fields_meta) &
                                  bind(c, name='fv3jedi_geom_addfmd_f90')

!Arguments
integer(c_int),     intent(inout) :: c_key_self
type(c_ptr), value, intent(in)    :: c_fields_meta

type(fv3jedi_geom), pointer :: self
type(fields_metadata)       :: f_fields_metadata

! LinkedList
! ----------
call fv3jedi_geom_registry%get(c_key_self,self)

! Fortran APIs
! ------------
f_fields_metadata = fields_metadata(c_fields_meta)

! Add to Fortran type
! -------------------
self%fmd = f_fields_metadata

end subroutine c_fv3jedi_addfmd_setup

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_clone(c_key_self, c_key_other, c_fields_meta) bind(c,name='fv3jedi_geom_clone_f90')

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

subroutine c_fv3jedi_geom_set_lonlat(c_key_self, c_afieldset, c_include_halo) &
                                     bind(c,name='fv3jedi_geom_set_lonlat_f90')

!Arguments
integer(c_int), intent(in) :: c_key_self
type(c_ptr), intent(in), value :: c_afieldset
logical(c_bool), intent(in) :: c_include_halo

type(fv3jedi_geom), pointer :: self
type(atlas_fieldset) :: afieldset
logical :: include_halo

! LinkedList
! ----------
call fv3jedi_geom_registry%get(c_key_self,self)

! Fortran APIs
! ------------
afieldset = atlas_fieldset(c_afieldset)

include_halo = c_include_halo

! Call implementation
! -------------------
call self%set_lonlat(afieldset, include_halo)

end subroutine c_fv3jedi_geom_set_lonlat

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_set_functionspace_pointer(c_key_self,c_afunctionspace,c_afunctionspace_for_bump) &
    bind(c,name='fv3jedi_geom_set_functionspace_pointer_f90')

integer(c_int), intent(in)     :: c_key_self
type(c_ptr), intent(in), value :: c_afunctionspace
type(c_ptr), intent(in), value :: c_afunctionspace_for_bump

type(fv3jedi_geom),pointer :: self

! LinkedList
! ----------
call fv3jedi_geom_registry%get(c_key_self,self)

! Create function space
! ---------------------
self%afunctionspace = atlas_functionspace(c_afunctionspace)
self%afunctionspace_for_bump = atlas_functionspace(c_afunctionspace_for_bump)

end subroutine c_fv3jedi_geom_set_functionspace_pointer

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_set_and_fill_geometry_fields(c_key_self, c_afieldset) &
                                       bind(c,name='fv3jedi_geom_set_and_fill_geometry_fields_f90')

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
call self%set_and_fill_geometry_fields(afieldset)

end subroutine c_fv3jedi_geom_set_and_fill_geometry_fields

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_start_end(c_key_self, ist, iend, jst, jend, kst, kend, npz) &
                                    bind(c, name='fv3jedi_geom_start_end_f90')

integer(c_int), intent( in) :: c_key_self
integer(c_int), intent(out) :: ist, iend, jst, jend, kst, kend, npz

type(fv3jedi_geom), pointer :: self

integer(c_int) :: itd ! iterator dimension

! LinkedList
! ----------
call fv3jedi_geom_registry%get(c_key_self, self)

itd = self%iterator_dimension

ist  = self%isc
iend = self%iec
jst  = self%jsc
jend = self%jec
! 3D iterator starts from 0 for surface variables
if (3 == itd) then
   kst = 0
else
   kst = 1
end if
kend = self%npz
npz  = self%npz

end subroutine c_fv3jedi_geom_start_end

! ------------------------------------------------------------------------------
!> C++ interface to get dimension of the GeometryIterator
subroutine c_fv3jedi_geom_iterator_dimension_f90(c_key_self, itd) &
                                           bind(c, name='fv3jedi_geom_iterator_dimension_f90')
  integer(c_int), intent( in) :: c_key_self
  integer(c_int), intent(out) :: itd ! iterator dimension

  type(fv3jedi_geom), pointer :: self
  call fv3jedi_geom_registry%get(c_key_self, self)

  itd = self%iterator_dimension
end subroutine c_fv3jedi_geom_iterator_dimension_f90

!--------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_verticalCoord(c_key_self, vc, npz, psurf) &
                                    bind(c, name='fv3jedi_geom_verticalCoord_f90')

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

subroutine c_fv3jedi_geom_get_data(c_key_self, c_nlev, c_ak, c_bk, c_ptop) &
           bind(c,name='fv3jedi_geom_get_data_f90')

! Arguments
integer(c_int), intent(in)  :: c_key_self
integer(c_int), intent(in)  :: c_nlev
real(c_double), intent(out) :: c_ak(c_nlev+1)
real(c_double), intent(out) :: c_bk(c_nlev+1)
real(c_double), intent(out) :: c_ptop

! Locals
type(fv3jedi_geom), pointer :: self
real(kind=kind_real) :: f_ak(c_nlev+1)
real(kind=kind_real) :: f_bk(c_nlev+1)
real(kind=kind_real) :: f_ptop

! LinkedList
! ----------
call fv3jedi_geom_registry%get(c_key_self, self)

! Call implementation
! -------------------
call self%get_data(f_ak, f_bk, f_ptop)

! Precision changes
! -----------------
c_ak = real(f_ak, kind=c_double)
c_bk = real(f_bk, kind=c_double)
c_ptop = real(f_ptop, kind=c_double)

end subroutine c_fv3jedi_geom_get_data

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_get_num_nodes_and_elements(c_key_self, c_num_nodes, c_num_tris, c_num_quads) &
    bind(c, name='fv3jedi_geom_get_num_nodes_and_elements_f90')
  integer(c_int), intent( in) :: c_key_self
  integer(c_int), intent(out) :: c_num_nodes
  integer(c_int), intent(out) :: c_num_tris
  integer(c_int), intent(out) :: c_num_quads

  integer :: num_nodes
  integer :: num_tris
  integer :: num_quads

  type(fv3jedi_geom), pointer :: self
  call fv3jedi_geom_registry%get(c_key_self, self)
  call self%get_num_nodes_and_elements(num_nodes, num_tris, num_quads)

  c_num_nodes = num_nodes
  c_num_tris = num_tris
  c_num_quads = num_quads

end subroutine c_fv3jedi_geom_get_num_nodes_and_elements

! --------------------------------------------------------------------------------------------------

subroutine c_fv3jedi_geom_get_coords_and_connectivities(c_key_self, &
    c_num_nodes, c_lons, c_lats, c_ghosts, c_global_indices, c_remote_indices, c_partition, &
    c_num_tri_boundary_nodes, c_raw_tri_boundary_nodes, &
    c_num_quad_boundary_nodes, c_raw_quad_boundary_nodes) &
    bind(c, name='fv3jedi_geom_get_coords_and_connectivities_f90')
  integer(c_int), intent( in) :: c_key_self
  integer(c_int), intent( in) :: c_num_nodes
  real(c_double), intent(out) :: c_lons(c_num_nodes)
  real(c_double), intent(out) :: c_lats(c_num_nodes)
  integer(c_int), intent(out) :: c_ghosts(c_num_nodes)
  integer(c_int), intent(out) :: c_global_indices(c_num_nodes)
  integer(c_int), intent(out) :: c_remote_indices(c_num_nodes)
  integer(c_int), intent(out) :: c_partition(c_num_nodes)
  integer(c_int), intent( in) :: c_num_tri_boundary_nodes
  integer(c_int), intent(out) :: c_raw_tri_boundary_nodes(c_num_tri_boundary_nodes)
  integer(c_int), intent( in) :: c_num_quad_boundary_nodes
  integer(c_int), intent(out) :: c_raw_quad_boundary_nodes(c_num_quad_boundary_nodes)

  integer :: num_nodes
  integer :: num_tri_boundary_nodes
  integer :: num_quad_boundary_nodes
  type(fv3jedi_geom), pointer :: self

  real(kind_real) :: lons(c_num_nodes)
  real(kind_real) :: lats(c_num_nodes)
  integer :: ghosts(c_num_nodes)
  integer :: global_indices(c_num_nodes)
  integer :: remote_indices(c_num_nodes)
  integer :: partition(c_num_nodes)
  integer :: raw_tri_boundary_nodes(c_num_tri_boundary_nodes)
  integer :: raw_quad_boundary_nodes(c_num_quad_boundary_nodes)

  num_nodes = c_num_nodes
  num_tri_boundary_nodes = c_num_tri_boundary_nodes
  num_quad_boundary_nodes = c_num_quad_boundary_nodes

  call fv3jedi_geom_registry%get(c_key_self, self)
  call self%get_coords_and_connectivities(num_nodes, num_tri_boundary_nodes, num_quad_boundary_nodes, &
    lons, lats, ghosts, global_indices, remote_indices, partition, raw_tri_boundary_nodes, raw_quad_boundary_nodes)

  c_lons = lons
  c_lats = lats
  c_ghosts = ghosts
  c_global_indices = global_indices
  c_remote_indices = remote_indices
  c_partition = partition
  c_raw_tri_boundary_nodes = raw_tri_boundary_nodes
  c_raw_quad_boundary_nodes = raw_quad_boundary_nodes

end subroutine c_fv3jedi_geom_get_coords_and_connectivities

! --------------------------------------------------------------------------------------------------

end module fv3jedi_geom_interface_mod
