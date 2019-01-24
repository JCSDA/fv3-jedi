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
  integer :: npx, npy, npz
  integer :: ntiles, ntile
  logical :: hydrostatic = .true.
  logical :: tladphystrj = .false.
  integer :: calendar_type, date_init(6) !Read/write for GFS
  integer :: nf

  type(fv3jedi_field), allocatable :: fields(:)

  integer :: ua   = 0
  integer :: va   = 0
  integer :: t    = 0
  integer :: ps   = 0
  integer :: q    = 0
  integer :: qi   = 0
  integer :: ql   = 0
  integer :: o3   = 0
  integer :: psi  = 0
  integer :: chi  = 0
  integer :: tv   = 0
  integer :: rh   = 0
  integer :: w    = 0
  integer :: delz = 0
  integer :: delp = 0

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
