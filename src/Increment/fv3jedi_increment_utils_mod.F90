! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Utilities for increment for the FV3JEDI model

module fv3jedi_increment_utils_mod

use kinds
use fv3jedi_vars_mod, only: fv3jedi_vars 

implicit none
private
public fv3jedi_increment

!> Fortran derived type to hold FV3JEDI increment
type :: fv3jedi_increment

  type(fv3jedi_vars) :: vars 

  !Local copies of grid
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: npx, npy, npz
  logical :: hydrostatic = .false.
  integer :: calendar_type
  integer, dimension(6) :: date
  integer, dimension(6) :: date_init

  !Note: for simplicity in transforming variables the increment is A-Grid winds.
  !This means all variables are co-located at cell centers.

  !Increment variables
  real(kind=kind_real), allocatable, dimension(:,:,:) :: ua     ! A-grid zonal wind (m/s)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: va     ! A-grid meridional wind (m/s)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: t      ! dry temperature (K)
  real(kind=kind_real), allocatable, dimension(:,:  ) :: ps     ! surface pressure (pascal)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: q      ! specific humidity (kg/kg)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: qi     ! cloud liquid ice (kg/kg)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: ql     ! cloud liquid water (kg/kg)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: o3     ! ozone (kg/kg)

  !Nonhydrostatic increment (not using for now)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: w      ! cell center vertical wind (m/s)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: delz   ! layer thickness (meters)

  !Control variables (used for the B-Matrix/Jb)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: psi    ! Stream function
  real(kind=kind_real), allocatable, dimension(:,:,:) :: chi    ! Velocity potential
  real(kind=kind_real), allocatable, dimension(:,:,:) :: tv     ! Virtual temperature
  real(kind=kind_real), allocatable, dimension(:,:,:) :: qc     ! humidity control variable 
  real(kind=kind_real), allocatable, dimension(:,:,:) :: qic    ! cloud liquid ice control variable
  real(kind=kind_real), allocatable, dimension(:,:,:) :: qlc    ! cloud liquid water control variable
  real(kind=kind_real), allocatable, dimension(:,:,:) :: o3c    ! ozone control variable

end type fv3jedi_increment

end module fv3jedi_increment_utils_mod
