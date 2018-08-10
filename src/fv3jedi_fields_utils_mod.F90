! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Utilities for fields for the FV3JEDI model

module fv3jedi_fields_utils_mod

use kinds
use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_vars_mod, only: fv3jedi_vars 

implicit none
private
public fv_atmos_type, fv3jedi_field

!> Skinny version of fv_atmos_type
type fv_atmos_type
  logical :: hydrostatic = .false.
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: u      ! A or D grid zonal wind (m/s)
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: v      ! A or D grid meridional wind (m/s)
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: pt     ! temperature (K)
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: delp   ! pressure thickness (pascal)
  real(kind=kind_real), allocatable, dimension(:,:,:,:) :: q      ! tracers (specific humidity and prognostic constituents)
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: w      ! cell center vertical wind (m/s)
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: delz   ! layer thickness (meters)
  real(kind=kind_real), allocatable, dimension(:,:)     :: phis   ! Surface geopotential (g*Z_surf)
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: ua     ! A grid zonal wind (m/s)
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: va     ! A grid meridional wind (m/s)

  !Control variables
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: vort   ! Vorticity
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: divg   ! Divergence
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: psi    ! Stream function
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: chi    ! Velocity potential
  real(kind=kind_real), allocatable, dimension(:,:,:)   :: tv     ! Virtual temperature
  real(kind=kind_real), allocatable, dimension(:,:)     :: ps     ! Surface temperature
  real(kind=kind_real), allocatable, dimension(:,:,:,:) :: qct    ! Tracer control variables

  integer :: calendar_type
  integer, dimension(6) :: date
  integer, dimension(6) :: date_init

  !2D Fields to be read
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
end type fv_atmos_type

!> Fortran derived type to hold FV3JEDI fields
type :: fv3jedi_field
  type(fv_atmos_type) :: Atm
  type(fv3jedi_vars) :: vars 
  type(fv3jedi_geom), pointer :: geom
!  integer :: nf
  integer :: isc, iec, jsc, jec, npz, nq
  integer :: isd, ied, jsd, jed
  integer :: root_pe
  logical :: havecrtmfields = .false.
  integer :: ti_q, ti_ql, ti_qi, ti_o3
end type fv3jedi_field

end module fv3jedi_fields_utils_mod
