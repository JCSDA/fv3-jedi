! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Utilities for state for the FV3JEDI model

module fv3jedi_state_utils_mod

use fv3jedi_kinds_mod
use fv3jedi_vars_mod, only: fv3jedi_vars 
use fv3jedi_field_mod, only: fv3jedi_field

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
  logical :: tladphystrj = .false.
  integer :: calendar_type, date_init(6) !Read/write for GFS
  integer :: nf

  type(fv3jedi_field), allocatable :: fields(:)

  !State variables (index in array for later access)
  integer :: ud   = 0
  integer :: vd   = 0
  integer :: ua   = 0
  integer :: va   = 0
  integer :: t    = 0
  integer :: delp = 0
  integer :: q    = 0
  integer :: qi   = 0
  integer :: ql   = 0
  integer :: o3   = 0
  integer :: w    = 0
  integer :: delz = 0
  integer :: phis = 0

  !CRTM state
  integer :: slmsk  = 0
  integer :: sheleg = 0
  integer :: tsea   = 0
  integer :: vtype  = 0
  integer :: stype  = 0
  integer :: vfrac  = 0
  integer :: stc    = 0
  integer :: smc    = 0
  integer :: snwdph = 0
  integer :: u_srf  = 0
  integer :: v_srf  = 0
  integer :: f10m   = 0

  !Linearized model trajectory
  integer :: qls     = 0
  integer :: qcn     = 0
  integer :: cfcn    = 0
  integer :: frocean = 0
  integer :: frland  = 0
  integer :: varflt  = 0
  integer :: ustar   = 0
  integer :: bstar   = 0
  integer :: zpbl    = 0
  integer :: cm      = 0
  integer :: ct      = 0
  integer :: cq      = 0
  integer :: kcbl    = 0
  integer :: ts      = 0
  integer :: khl     = 0
  integer :: khu     = 0

end type fv3jedi_state

! ------------------------------------------------------------------------------

#define LISTED_TYPE fv3jedi_state

!> Linked list interface - defines registry_t type
#include "linkedList_i.f"

!> Global registry
type(registry_t) :: fv3jedi_state_registry

! ------------------------------------------------------------------------------

contains

! ------------------------------------------------------------------------------

!> Linked list implementation
#include "linkedList_c.f"

! ------------------------------------------------------------------------------

end module fv3jedi_state_utils_mod
