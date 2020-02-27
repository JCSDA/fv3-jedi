! (C) Copyright 2018 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> Utilities for state for the FV3JEDI model

module fv3jedi_state_utils_mod

use fv3jedi_kinds_mod
use fv3jedi_field_mod, only: fv3jedi_field
use fckit_mpi_module, only: fckit_mpi_comm

implicit none
private
public fv3jedi_state, fv3jedi_state_registry

!> Fortran derived type to hold FV3JEDI state
type :: fv3jedi_state

  !Local copies of grid for convenience
  integer :: isc, iec, jsc, jec
  integer :: npx, npy, npz
  integer :: ntiles, ntile
  logical :: hydrostatic = .true.
  integer :: calendar_type, date_init(6) !Read/write for GFS
  integer :: nf
  logical :: have_agrid
  logical :: have_dgrid

  type(fckit_mpi_comm) :: f_comm
  type(fv3jedi_field), allocatable :: fields(:)

end type fv3jedi_state

! ------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_state

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_state_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

end module fv3jedi_state_utils_mod
