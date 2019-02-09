! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Utilities for increment for the FV3JEDI model

module fv3jedi_increment_utils_mod

use fv3jedi_kinds_mod
use fv3jedi_field_mod, only: fv3jedi_field

implicit none
private
public fv3jedi_increment, fv3jedi_increment_registry

!> Fortran derived type to hold FV3JEDI increment
type :: fv3jedi_increment

  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: npx, npy, npz
  integer :: ntiles, ntile
  logical :: hydrostatic = .true.
  logical :: tladphystrj = .false.
  integer :: calendar_type, date_init(6) !Read/write for GFS
  integer :: nf

  type(fv3jedi_field), allocatable :: fields(:)

  real(kind=kind_real), pointer, dimension(:,:,:) :: ua   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: va   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: t    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ps   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: q    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: qi   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ql   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: o3   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: psi  => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: chi  => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: tv   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: rh   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: w    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: delz => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: delp => null()

end type fv3jedi_increment

! ------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_increment

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_increment_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

end module fv3jedi_increment_utils_mod
