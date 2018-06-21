! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Utilities for fields for the FV3JEDI model

module fv3jedi_fields_utils_mod

use kinds
use mpp_mod,         only: mpp_pe, mpp_npes, mpp_error, FATAL, NOTE
use mpp_domains_mod, only: domain2D, mpp_define_layout, mpp_define_mosaic, mpp_define_io_domain
use fv3jedi_geom_mod,only: fv3jedi_geom

implicit none
private
public fv_atmos_type, fv3jedi_field, allocate_fv_atmos_type

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
  type(fv3jedi_geom), pointer :: geom
  integer :: nf
  integer :: isc, iec, jsc, jec
  integer :: root_pe
  logical :: havecrtmfields
end type fv3jedi_field

contains

! Allocate the main model/increment fields
! ----------------------------------------
subroutine allocate_fv_atmos_type(Atm, isd, ied, jsd, jed, &
                                       isc, iec, jsc, jec, &
                                       nz, nq, hydrostatic, wind_type)

 implicit none

 type(fv_atmos_type), intent(inout), target :: Atm
 logical, intent(in) :: hydrostatic
 integer, intent(in) :: isd, ied, jsd, jed
 integer, intent(in) :: isc, iec, jsc, jec
 integer, intent(in) :: nz, nq
 character(len=255)  :: wind_type

  Atm%hydrostatic = hydrostatic

  if (trim(wind_type) == 'D-grid') then
     if (.not.allocated(   Atm%u)) allocate (    Atm%u(isd:ied,   jsd:jed+1 , nz     ) )
     if (.not.allocated(   Atm%v)) allocate (    Atm%v(isd:ied+1, jsd:jed   , nz     ) )
  elseif (trim(wind_type) == 'A-grid') then
     if (.not.allocated(   Atm%u)) allocate (    Atm%u(isd:ied, jsd:jed, nz     ) )
     if (.not.allocated(   Atm%v)) allocate (    Atm%v(isd:ied, jsd:jed, nz     ) )
  else
     call abor1_ftn("module_fv3jedi: wind_type issue, must be either A-grid or D-grid")
  endif

  if (.not.allocated(  Atm%pt)) allocate (   Atm%pt(isd:ied,   jsd:jed   , nz     ) )
  if (.not.allocated(Atm%delp)) allocate ( Atm%delp(isd:ied,   jsd:jed   , nz     ) )
  if (.not.allocated(   Atm%q)) allocate (    Atm%q(isd:ied,   jsd:jed   , nz, nq ) )
  if (.not.allocated(Atm%phis)) allocate ( Atm%phis(isd:ied,   jsd:jed            ) )

  if (.not. hydrostatic) then
     if (.not.allocated(   Atm%w)) allocate (    Atm%w(isd:ied,   jsd:jed   , nz     ) )
     if (.not.allocated(Atm%delz)) allocate ( Atm%delz(isd:ied,   jsd:jed   , nz     ) )
  endif

  if (.not.allocated(   Atm%ua)) allocate (    Atm%ua(isd:ied, jsd:jed, nz     ) )
  if (.not.allocated(   Atm%va)) allocate (    Atm%va(isd:ied, jsd:jed, nz     ) )

  if (.not.allocated(Atm%slmsk )) allocate(Atm%slmsk (isd:ied,jsd:jed))
  if (.not.allocated(Atm%sheleg)) allocate(Atm%sheleg(isd:ied,jsd:jed))
  if (.not.allocated(Atm%tsea  )) allocate(Atm%tsea  (isd:ied,jsd:jed))
  if (.not.allocated(Atm%vtype )) allocate(Atm%vtype (isd:ied,jsd:jed))
  if (.not.allocated(Atm%stype )) allocate(Atm%stype (isd:ied,jsd:jed))
  if (.not.allocated(Atm%vfrac )) allocate(Atm%vfrac (isd:ied,jsd:jed))
  if (.not.allocated(Atm%stc   )) allocate(Atm%stc   (isd:ied,jsd:jed,4))
  if (.not.allocated(Atm%smc   )) allocate(Atm%smc   (isd:ied,jsd:jed,4))
  if (.not.allocated(Atm%u_srf )) allocate(Atm%snwdph(isd:ied,jsd:jed))
  if (.not.allocated(Atm%u_srf )) allocate(Atm%u_srf (isd:ied,jsd:jed))
  if (.not.allocated(Atm%v_srf )) allocate(Atm%v_srf (isd:ied,jsd:jed))
  if (.not.allocated(Atm%f10m  )) allocate(Atm%f10m  (isd:ied,jsd:jed))

end subroutine allocate_fv_atmos_type

end module fv3jedi_fields_utils_mod
