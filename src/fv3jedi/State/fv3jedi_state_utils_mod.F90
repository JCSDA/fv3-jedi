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

  !Adding a new variable? Update create, delete and copy
  real(kind=kind_real), pointer, dimension(:,:,:) :: ud      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: vd      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ua      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: va      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: t       => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: tv      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: pt      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: delp    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: pe      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: pkz     => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ps      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: q       => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: rh      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: qi      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ql      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: qils    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: qlls    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: qicn    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: qlcn    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: qs      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: qr      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: gr      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ca      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: o3      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ox      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: w       => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: delz    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: phis    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: psi     => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: chi     => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: vort    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: divg    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: slmsk   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: sheleg  => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: tsea    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: vtype   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: stype   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: vfrac   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: stc     => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: smc     => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: snwdph  => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: u_srf   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: v_srf   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: f10m    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: qls     => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: qcn     => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: cfcn    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: frocean => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: frland  => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: varflt  => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ustar   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: bstar   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: zpbl    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: cm      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ct      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: cq      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: kcbl    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ts      => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: khl     => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: khu     => null()
  !Aerosols
  real(kind=kind_real), pointer, dimension(:,:,:) :: du001    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: du002    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: du003    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: du004    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: du005    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ss001    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ss002    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ss003    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ss004    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ss005    => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: no3an1   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: no3an2   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: no3an3   => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: bcphobic => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: bcphilic => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ocphobic => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: ocphilic => null()
  real(kind=kind_real), pointer, dimension(:,:,:) :: so4      => null()

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
