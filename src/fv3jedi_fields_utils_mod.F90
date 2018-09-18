! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Utilities for fields for the FV3JEDI model

module fv3jedi_fields_utils_mod

use kinds
use fv3jedi_vars_mod, only: fv3jedi_vars 

implicit none
private
public fv3jedi_field

!> Fortran derived type to hold FV3JEDI fields
type :: fv3jedi_field

  type(fv3jedi_vars) :: vars 

  !Local copies of grid
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: npx, npy, npz
  logical :: havecrtmfields = .false.
  logical :: hydrostatic = .false.
  integer :: calendar_type
  integer, dimension(6) :: date
  integer, dimension(6) :: date_init
  logical :: am_i_root_pe 

  !State/increment variables
  real(kind=kind_real), allocatable, dimension(:,:,:) :: ud     ! D-grid (grid tangential) zonal wind (m/s)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: vd     ! D-grid (grid tangential) meridional wind (m/s)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: ua     ! A-grid zonal wind (m/s)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: va     ! A-grid meridional wind (m/s)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: t      ! dry temperature (K)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: delp   ! pressure thickness (pascal)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: q      ! specific humidity (kg/kg)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: qi     ! cloud liquid ice (kg/kg)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: ql     ! cloud liquid water (kg/kg)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: o3     ! ozone (kg/kg)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: w      ! cell center vertical wind (m/s)
  real(kind=kind_real), allocatable, dimension(:,:,:) :: delz   ! layer thickness (meters)
  real(kind=kind_real), allocatable, dimension(:,:)   :: phis   ! Surface geopotential (g*Z_surf)

  !Control/analysis variables
  real(kind=kind_real), allocatable, dimension(:,:,:) :: vort   ! Vorticity
  real(kind=kind_real), allocatable, dimension(:,:,:) :: divg   ! Divergence
  real(kind=kind_real), allocatable, dimension(:,:,:) :: psi    ! Stream function
  real(kind=kind_real), allocatable, dimension(:,:,:) :: chi    ! Velocity potential
  real(kind=kind_real), allocatable, dimension(:,:,:) :: tv     ! Virtual temperature
  real(kind=kind_real), allocatable, dimension(:,:)   :: ps     ! Surface temperature
  real(kind=kind_real), allocatable, dimension(:,:,:) :: qc     ! humidity control variable 
  real(kind=kind_real), allocatable, dimension(:,:,:) :: qic    ! cloud liquid ice control variable
  real(kind=kind_real), allocatable, dimension(:,:,:) :: qlc    ! cloud liquid water control variable
  real(kind=kind_real), allocatable, dimension(:,:,:) :: o3c    ! ozone control variable

  !2D Fields for CRTM
  integer             , allocatable, dimension(:,:)   :: slmsk
  real(kind=kind_real), allocatable, dimension(:,:)   :: sheleg
  real(kind=kind_real), allocatable, dimension(:,:)   :: tsea
  integer             , allocatable, dimension(:,:)   :: vtype
  integer             , allocatable, dimension(:,:)   :: stype
  real(kind=kind_real), allocatable, dimension(:,:)   :: vfrac
  real(kind=kind_real), allocatable, dimension(:,:,:) :: stc
  real(kind=kind_real), allocatable, dimension(:,:,:) :: smc
  real(kind=kind_real), allocatable, dimension(:,:)   :: snwdph
  real(kind=kind_real), allocatable, dimension(:,:)   :: u_srf
  real(kind=kind_real), allocatable, dimension(:,:)   :: v_srf
  real(kind=kind_real), allocatable, dimension(:,:)   :: f10m

end type fv3jedi_field

end module fv3jedi_fields_utils_mod
