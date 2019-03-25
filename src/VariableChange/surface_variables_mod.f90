! (C) Copyright 2018-2019 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Variable transforms/interpolation for surface variables in fv3-jedi 
!> Daniel Holdaway, NASA/JCSDA

module surface_vt_mod

use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_kinds_mod, only: kind_real
use crtm_module, only: crtm_irlandcoeff_classification
use fv3jedi_constants_mod, only: rad2deg, deg2rad, pi
use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum
use type_kdtree, only: kdtree_type
use type_mpl, only: mpl_type
use tools_func, only: sphere_dist

implicit none
private

public crtm_surface

contains

!----------------------------------------------------------------------------
! Surface quantities in the form needed by the crtm ------------------------- 
!----------------------------------------------------------------------------

subroutine crtm_surface( geom, nobs, ngrid, lats_ob, lons_ob, &
                         fld_slmsk, fld_sheleg, fld_tsea, fld_vtype, fld_stype, fld_vfrac, fld_stc, &
                         fld_smc, fld_snwdph, fld_u_srf, fld_v_srf, fld_f10m, &
                         land_type, vegetation_type, soil_type, water_coverage, land_coverage, ice_coverage, &
                         snow_coverage, lai, water_temperature, land_temperature, ice_temperature, &
                         snow_temperature, soil_moisture_content, vegetation_fraction, soil_temperature, snow_depth, &
                         wind_speed, wind_direction )

implicit none

!Arguments
type(fv3jedi_geom)  , intent(in)  :: geom !Geometry for the model
integer             , intent(in)  :: nobs
integer             , intent(in)  :: ngrid
real(kind=kind_real), intent(in)  :: lats_ob(nobs)
real(kind=kind_real), intent(in)  :: lons_ob(nobs)
real(kind=kind_real), intent(in)  :: fld_slmsk (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: fld_sheleg(geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: fld_tsea  (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: fld_vtype (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: fld_stype (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: fld_vfrac (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: fld_stc   (geom%isc:geom%iec,geom%jsc:geom%jec,4)
real(kind=kind_real), intent(in)  :: fld_smc   (geom%isc:geom%iec,geom%jsc:geom%jec,4)
real(kind=kind_real), intent(in)  :: fld_snwdph(geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: fld_u_srf (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: fld_v_srf (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: fld_f10m  (geom%isc:geom%iec,geom%jsc:geom%jec,1)
integer             , intent(out) :: vegetation_type(nobs)
integer             , intent(out) :: land_type(nobs)
integer             , intent(out) :: soil_type(nobs)
real(kind=kind_real), intent(out) :: water_coverage(nobs)
real(kind=kind_real), intent(out) :: land_coverage(nobs)
real(kind=kind_real), intent(out) :: ice_coverage(nobs)
real(kind=kind_real), intent(out) :: snow_coverage(nobs)
real(kind=kind_real), intent(out) :: lai(nobs)
real(kind=kind_real), intent(out) :: water_temperature(nobs)
real(kind=kind_real), intent(out) :: land_temperature(nobs)
real(kind=kind_real), intent(out) :: ice_temperature(nobs)
real(kind=kind_real), intent(out) :: snow_temperature(nobs)
real(kind=kind_real), intent(out) :: soil_moisture_content(nobs)
real(kind=kind_real), intent(out) :: vegetation_fraction(nobs)
real(kind=kind_real), intent(out) :: soil_temperature(nobs)
real(kind=kind_real), intent(out) :: snow_depth(nobs)
real(kind=kind_real), intent(out) :: wind_speed(nobs)
real(kind=kind_real), intent(out) :: wind_direction(nobs)

!Locals
real(kind=kind_real), parameter :: minsnow = 1.0_kind_real / 10.0_kind_real
real(kind=kind_real), parameter :: windlimit = 0.0001_kind_real
real(kind=kind_real), parameter :: quadcof  (4, 2  ) =      &
                                   reshape((/0.0_kind_real,  1.0_kind_real, 1.0_kind_real,  2.0_kind_real, &
                                             1.0_kind_real, -1.0_kind_real, 1.0_kind_real, -1.0_kind_real/), (/4, 2/))

integer              :: n, itype, istype
integer              :: istyp00,istyp01,istyp10,istyp11
integer              :: isflg, idomsfc
integer              :: lai_type, iquadrant
logical              :: lwind
real(kind=kind_real) :: dtsfc, dtsfcp
real(kind=kind_real) :: sfcpct(0:3), ts(0:3), wgtavg(0:3), dtskin(0:3)
real(kind=kind_real) :: w00,w01,w10,w11
real(kind=kind_real) :: sno00,sno01,sno10,sno11
real(kind=kind_real) :: sst00,sst01,sst10,sst11
real(kind=kind_real) :: tsavg,wgtmin
real(kind=kind_real) :: vty, sty, vfr, stp, sm, sn
real(kind=kind_real) :: uu5, vv5, f10, sfc_speed, windratio, windangle, windscale
real(kind=kind_real) :: wind10, wind10_direction

!From GSI
integer, parameter :: GFS_SOIL_N_TYPES = 9
integer, parameter :: GFS_VEGETATION_N_TYPES = 13
integer, parameter :: INVALID_LAND = 0
integer, parameter :: COMPACTED_SOIL = 1
integer, parameter :: TILLED_SOIL = 2
integer, parameter :: IRRIGATED_LOW_VEGETATION = 5
integer, parameter :: MEADOW_GRASS = 6
integer, parameter :: SCRUB = 7
integer, parameter :: BROADLEAF_FOREST = 8
integer, parameter :: PINE_FOREST = 9
integer, parameter :: TUNDRA = 10
integer, parameter :: GRASS_SOIL = 11
integer, parameter :: BROADLEAF_PINE_FOREST = 12
integer, parameter :: GRASS_SCRUB = 13
integer, parameter :: URBAN_CONCRETE = 15
integer, parameter :: BROADLEAF_BRUSH = 17
integer, parameter :: WET_SOIL = 18
integer, parameter :: SCRUB_SOIL = 19
integer, parameter :: nvege_type = 20
integer, parameter :: IGBP_N_TYPES = 20
integer, parameter :: SOIL_N_TYPES = 16
integer, allocatable,dimension(:) :: map_to_crtm_ir
integer, allocatable,dimension(:) :: map_to_crtm_mwave
integer, parameter, dimension(1:IGBP_N_TYPES) :: igbp_to_gfs=(/4, &
   1, 5, 2, 3, 8, 9, 6, 6, 7, 8, 12, 7, 12, 13, 11, 0, 10, 10, 11/)
 integer, parameter, dimension(1:IGBP_N_TYPES) :: igbp_to_npoess=(/PINE_FOREST, &
   BROADLEAF_FOREST, PINE_FOREST, BROADLEAF_FOREST, BROADLEAF_PINE_FOREST, &
   SCRUB, SCRUB_SOIL, BROADLEAF_BRUSH, BROADLEAF_BRUSH, SCRUB, BROADLEAF_BRUSH, &
   TILLED_SOIL, URBAN_CONCRETE, TILLED_SOIL, INVALID_LAND, COMPACTED_SOIL, &
   INVALID_LAND, TUNDRA, TUNDRA, TUNDRA/)
 integer, parameter, dimension(1:IGBP_N_TYPES) :: igbp_to_igbp=(/1, &
   2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, &
   20/)
 integer, parameter, dimension(1:SOIL_N_TYPES) :: map_soil_to_crtm=(/1, &
   1, 4, 2, 2, 8, 7, 2, 6, 5, 2, 3, 8, 1, 6, 9/)

!For interpolation
integer :: nn
real(kind=kind_real), allocatable, dimension(:,:) :: interp_w
integer             , allocatable, dimension(:,:) :: interp_i

integer             , allocatable, dimension(:,:) :: slmsk
real(kind=kind_real), allocatable, dimension(:,:) :: sheleg
real(kind=kind_real), allocatable, dimension(:,:) :: tsea
integer             , allocatable, dimension(:,:) :: vtype
integer             , allocatable, dimension(:,:) :: stype
real(kind=kind_real), allocatable, dimension(:,:) :: vfrac
real(kind=kind_real), allocatable, dimension(:,:) :: snwdph
real(kind=kind_real), allocatable, dimension(:,:) :: stc
real(kind=kind_real), allocatable, dimension(:,:) :: smc
real(kind=kind_real), allocatable, dimension(:,:) :: u_srf
real(kind=kind_real), allocatable, dimension(:,:) :: v_srf
real(kind=kind_real), allocatable, dimension(:,:) :: f10m

real(kind=kind_real), allocatable, dimension(:,:) :: rslmsk
real(kind=kind_real), allocatable, dimension(:,:) :: rvtype
real(kind=kind_real), allocatable, dimension(:,:) :: rstype

integer             , allocatable, dimension(:,:) :: slmskp
real(kind=kind_real), allocatable, dimension(:,:) :: shelegp
real(kind=kind_real), allocatable, dimension(:,:) :: tseap
integer             , allocatable, dimension(:,:) :: vtypep
integer             , allocatable, dimension(:,:) :: stypep
real(kind=kind_real), allocatable, dimension(:,:) :: vfracp
real(kind=kind_real), allocatable, dimension(:,:) :: snwdphp
real(kind=kind_real), allocatable, dimension(:,:) :: stcp
real(kind=kind_real), allocatable, dimension(:,:) :: smcp
real(kind=kind_real), allocatable, dimension(:,:) :: u_srfp
real(kind=kind_real), allocatable, dimension(:,:) :: v_srfp
real(kind=kind_real), allocatable, dimension(:,:) :: f10mp

type(fckit_mpi_comm) :: comm

! Communicator
comm = fckit_mpi_comm() 

! Number of weights
nn = 4

allocate(interp_w(nobs,nn))
allocate(interp_i(nobs,nn))

allocate(slmsk(nobs,nn))
allocate(sheleg(nobs,nn))
allocate(tsea(nobs,nn))
allocate(vtype(nobs,nn))
allocate(stype(nobs,nn))
allocate(vfrac(nobs,nn))
allocate(stc(nobs,nn))
allocate(smc(nobs,nn))
allocate(snwdph(nobs,nn))
allocate(u_srf(nobs,nn))
allocate(v_srf(nobs,nn))
allocate(f10m(nobs,nn))

allocate(rslmsk(nobs,nn))
allocate(rvtype(nobs,nn))
allocate(rstype(nobs,nn))

allocate(slmskp(nobs,nn))
allocate(shelegp(nobs,nn))
allocate(tseap(nobs,nn))
allocate(vtypep(nobs,nn))
allocate(stypep(nobs,nn))
allocate(vfracp(nobs,nn))
allocate(stcp(nobs,nn))
allocate(smcp(nobs,nn))
allocate(snwdphp(nobs,nn))
allocate(u_srfp(nobs,nn))
allocate(v_srfp(nobs,nn))
allocate(f10mp(nobs,nn))

 !Second time level option, zero for now
 slmskp = 0
 shelegp = 0.0_kind_real
 tseap = 0.0_kind_real
 vtypep = 0.0_kind_real
 stypep = 0
 vfracp = 0
 snwdphp = 0.0_kind_real
 stcp = 0.0_kind_real
 smcp = 0.0_kind_real
 snwdphp = 0.0_kind_real
 u_srfp = 0.0_kind_real
 v_srfp = 0.0_kind_real
 f10mp = 0.0_kind_real

 !Get interpolation weights and index in global grid
 call kdtree_getweightsindex(comm, geom, nobs, ngrid, deg2rad*lats_ob, deg2rad*lons_ob, &
                             nn, interp_w, interp_i) 

 !Get field at nn neighbours
 call kdtree_getneighbours(comm, geom, nobs, ngrid, nn, interp_i, fld_slmsk (:,:,1), rslmsk)
 call kdtree_getneighbours(comm, geom, nobs, ngrid, nn, interp_i, fld_sheleg(:,:,1), sheleg)
 call kdtree_getneighbours(comm, geom, nobs, ngrid, nn, interp_i, fld_tsea  (:,:,1), tsea  )
 call kdtree_getneighbours(comm, geom, nobs, ngrid, nn, interp_i, fld_vtype (:,:,1), rvtype)
 call kdtree_getneighbours(comm, geom, nobs, ngrid, nn, interp_i, fld_stype (:,:,1), rstype)
 call kdtree_getneighbours(comm, geom, nobs, ngrid, nn, interp_i, fld_vfrac (:,:,1), vfrac )
 call kdtree_getneighbours(comm, geom, nobs, ngrid, nn, interp_i, fld_stc   (:,:,1), stc   )
 call kdtree_getneighbours(comm, geom, nobs, ngrid, nn, interp_i, fld_smc   (:,:,1), smc   )
 call kdtree_getneighbours(comm, geom, nobs, ngrid, nn, interp_i, fld_snwdph(:,:,1), snwdph)
 call kdtree_getneighbours(comm, geom, nobs, ngrid, nn, interp_i, fld_u_srf (:,:,1), u_srf )
 call kdtree_getneighbours(comm, geom, nobs, ngrid, nn, interp_i, fld_v_srf (:,:,1), v_srf )
 call kdtree_getneighbours(comm, geom, nobs, ngrid, nn, interp_i, fld_f10m  (:,:,1), f10m  )

 !Convert some reals to integer
 slmsk = nint(rslmsk)
 stype = nint(rvtype)
 stype = nint(rstype)

 dtsfc  = 1.0_kind_real
 dtsfcp = 0.0_kind_real

 dtskin = 0.0_kind_real !TODO need real skin temperature increment?

 lwind = .true.

! Vegation, land and soil types
! -----------------------------
 allocate(map_to_crtm_ir(nvege_type))
 allocate(map_to_crtm_mwave(nvege_type))

 map_to_crtm_mwave=igbp_to_gfs
 map_to_crtm_ir = igbp_to_igbp

 !TODO, get this from CRTM properly
 !select case ( TRIM(CRTM_IRlandCoeff_Classification()) )
 ! case('NPOESS'); map_to_crtm_ir=igbp_to_npoess
 ! case('IGBP')  ; map_to_crtm_ir=igbp_to_igbp
 !end select

 !Loop over all obs
 
 do n = 1,nobs
 
! Stage 1, like deter_sfc in GSI
! ------------------------------

    w00 = interp_w(n,1)
    w10 = interp_w(n,2)
    w01 = interp_w(n,3)
    w11 = interp_w(n,4)
 
    istyp00 = slmsk(n,1)
    istyp10 = slmsk(n,2)
    istyp01 = slmsk(n,3)
    istyp11 = slmsk(n,4)
 
    sno00 = snwdph(n,1)*dtsfc + snwdphp(n,1)*dtsfcp
    sno01 = snwdph(n,2)*dtsfc + snwdphp(n,2)*dtsfcp
    sno10 = snwdph(n,3)*dtsfc + snwdphp(n,3)*dtsfcp
    sno11 = snwdph(n,4)*dtsfc + snwdphp(n,4)*dtsfcp

    sst00 = tsea(n,1)*dtsfc + tsea(n,1)*dtsfcp
    sst01 = tsea(n,2)*dtsfc + tsea(n,2)*dtsfcp
    sst10 = tsea(n,3)*dtsfc + tsea(n,3)*dtsfcp
    sst11 = tsea(n,4)*dtsfc + tsea(n,4)*dtsfcp
 
    tsavg = sst00*w00 + sst10*w10 + sst01*w01 + sst11*w11
 
    if (istyp00 >=1 .and. sno00 > minsnow) istyp00 = 3
    if (istyp01 >=1 .and. sno01 > minsnow) istyp01 = 3
    if (istyp10 >=1 .and. sno10 > minsnow) istyp10 = 3
    if (istyp11 >=1 .and. sno11 > minsnow) istyp11 = 3
 
    sfcpct = 0.0_kind_real
    sfcpct(istyp00) = sfcpct(istyp00) + w00
    sfcpct(istyp01) = sfcpct(istyp01) + w01
    sfcpct(istyp10) = sfcpct(istyp10) + w10
    sfcpct(istyp11) = sfcpct(istyp11) + w11
 
    isflg = 0
    if(sfcpct(0) > 0.99_kind_real)then
       isflg = 0
    else if(sfcpct(1) > 0.99_kind_real)then
       isflg = 1
    else if(sfcpct(2) > 0.99_kind_real)then
       isflg = 2
    else if(sfcpct(3) > 0.99_kind_real)then
       isflg = 3
    else
       isflg = 4
    end if
 
    ts(0:3)=0.0_kind_real
    wgtavg(0:3)=0.0_kind_real
    vfr=0.0_kind_real
    stp=0.0_kind_real
    sty=0.0_kind_real
    vty=0.0_kind_real
    sm=0.0_kind_real
    sn=0.0_kind_real
 
    idomsfc=slmsk(n,1)
    wgtmin = w00

    if(istyp00 == 1)then
       vty  = vtype(n,1)
       sty  = stype(n,1)
       wgtavg(1) = wgtavg(1) + w00
       ts(1)=ts(1)+w00*sst00
       vfr  =vfr  +w00*( vfrac(n,1) * dtsfc + &
                        vfracp(n,1) * dtsfcp  )
       stp  =stp  +w00*(   stc(n,1) * dtsfc + &
                          stcp(n,1) * dtsfcp  )
       sm   =sm   +w00*(   smc(n,1) * dtsfc + &
                          smcp(n,1) * dtsfcp  )
    else if(istyp00 == 2)then
       wgtavg(2) = wgtavg(2) + w00
       ts(2)=ts(2)+w00*sst00
    else if(istyp00 == 3)then
       wgtavg(3) = wgtavg(3) + w00
       ts(3)=ts(3)+w00*sst00
       sn = sn + w00*sno00
    else
       wgtavg(0) = wgtavg(0) + w00
       ts(0)=ts(0)+w00*sst00
    end if

    if(istyp01 == 1)then
       if(wgtmin < w01 .or. (vty == 0.0_kind_real .and. sty == 0.0_kind_real))then
          vty  = vtype(n,3)
          sty  = stype(n,3)
       end if
       wgtavg(1) = wgtavg(1) + w01
       ts(1)=ts(1)+w01*sst01
       vfr  =vfr  +w01*( vfrac(n,3) * dtsfc + &
                        vfracp(n,3) * dtsfcp  )
       stp  =stp  +w01*(   stc(n,3) * dtsfc + &
                          stcp(n,3) * dtsfcp  )
       sm   =sm   +w01*(   smc(n,3) * dtsfc + &
                          smcp(n,3) * dtsfcp  )
    else if(istyp01 == 2)then
       wgtavg(2) = wgtavg(2) + w01
       ts(2)=ts(2)+w01*sst01
    else if(istyp01 == 3)then
       wgtavg(3) = wgtavg(3) + w01
       ts(3)=ts(3)+w01*sst01
       sn = sn + w01*sno01
    else
       wgtavg(0) = wgtavg(0) + w01
       ts(0)=ts(0)+w01*sst01
    end if
    if(wgtmin < w01)then
       idomsfc=slmsk(n,3)
       wgtmin = w01
    end if

    if(istyp10 == 1)then
       if(wgtmin < w10 .or. (vty == 0.0_kind_real .and. sty == 0.0_kind_real))then
          vty  = vtype(n,2)
          sty  = stype(n,2)
       end if
       wgtavg(1) = wgtavg(1) + w10
       ts(1)=ts(1)+w10*sst10
       vfr  =vfr  +w10*(vfrac (n,2) * dtsfc + &
                        vfracp(n,2) * dtsfcp  )
       stp  =stp  +w10*(  stc (n,2) * dtsfc + &
                          stcp(n,2) * dtsfcp  )
       sm   =sm   +w10*(  smc (n,2) * dtsfc + &
                          smcp(n,2) * dtsfcp  )
    else if(istyp10 == 2)then
       wgtavg(2) = wgtavg(2) + w10
       ts(2)=ts(2)+w10*sst10
    else if(istyp10 == 3)then
       wgtavg(3) = wgtavg(3) + w10
       ts(3)=ts(3)+w10*sst10
       sn = sn + w10*sno10
    else
       wgtavg(0) = wgtavg(0) + w10
       ts(0)=ts(0)+w10*sst10
    end if
    if(wgtmin < w10)then
       idomsfc=slmsk(n,2)
       wgtmin = w10
    end if

    if(istyp11 == 1)then
       if(wgtmin < w11 .or. (vty == 0.0_kind_real .and. sty == 0.0_kind_real))then
          vty  = vtype(n,4)
          sty  = stype(n,4)
       endif
       wgtavg(1) = wgtavg(1) + w11
       ts(1)=ts(1)+w11*sst11
       vfr  =vfr  +w11*(vfrac (n,4) * dtsfc + &
                        vfracp(n,4) * dtsfcp  )
       stp  =stp  +w11*(  stc (n,4) * dtsfc + &
                          stcp(n,4) * dtsfcp  )
       sm   =sm   +w11*(  smc (n,4) * dtsfc + &
                          smcp(n,4) * dtsfcp  )
    else if(istyp11 == 2)then
       wgtavg(2) = wgtavg(2) + w11
       ts(2)=ts(2)+w11*sst11
    else if(istyp11 == 3)then
       wgtavg(3) = wgtavg(3) + w11
       ts(3)=ts(3)+w11*sst11
       sn = sn + w11*sno11
    else
       wgtavg(0) = wgtavg(0) + w11
       ts(0)=ts(0)+w11*sst11
    end if

    if(wgtmin < w11)then
       idomsfc=slmsk(n,4)
       wgtmin = w11
    end if
 
    if(wgtavg(0) > 0.0_kind_real)then
       ts(0) = ts(0)/wgtavg(0)
    else
       ts(0) = tsavg
    end if

    if(wgtavg(1) > 0.0_kind_real)then
       ts(1) = ts(1)/wgtavg(1)
       sm = sm/wgtavg(1)
       vfr = vfr/wgtavg(1)
       stp = stp/wgtavg(1)
    else
       ts(1) = tsavg
       sm=1.0_kind_real
    end if

    if(wgtavg(2) > 0.0_kind_real)then
       ts(2) = ts(2)/wgtavg(2)
    else
       ts(2) = tsavg
    end if

    if(wgtavg(3) > 0.0_kind_real)then
       ts(3) = ts(3)/wgtavg(3)
       sn = sn/wgtavg(3)
    else
       ts(3) = tsavg
    end if

    f10 = ( f10m (n,1)*w00 + f10m (n,2)*w10 + f10m (n,3)*w01 + f10m (n,4)*w11 ) * dtsfc + &
          ( f10mp(n,1)*w00 + f10mp(n,2)*w10 + f10mp(n,3)*w01 + f10mp(n,4)*w11 ) * dtsfcp

! Stage 2 - like crtm_interface from GSI
! --------------------------------------

   itype  = vty
   istype = sty

   itype  = min(max(1,itype),nvege_type)
   istype = min(max(1,istype),SOIL_N_TYPES)
   land_type(n) = max(1,map_to_crtm_ir(itype))
   Vegetation_Type(n) = max(1,map_to_crtm_mwave(itype))
   Soil_Type(n) = map_soil_to_crtm(istype)
   lai_type = map_to_crtm_mwave(itype)

   water_coverage(n) = min(max(0.0_kind_real,sfcpct(0)),1.0_kind_real)
   land_coverage(n)  = min(max(0.0_kind_real,sfcpct(1)),1.0_kind_real)
   ice_coverage(n)   = min(max(0.0_kind_real,sfcpct(2)),1.0_kind_real)
   snow_coverage(n)  = min(max(0.0_kind_real,sfcpct(3)),1.0_kind_real)

   Lai(n) = 0.0_kind_real

   if (land_coverage(n) > 0.0_kind_real) then

      if(lai_type>0)then
        call get_lai(lai_type,lai(n)) !TODO: does nothing yet
      endif     
   
      ! for Glacial land ice soil type and vegetation type
      if(Soil_Type(n) == 9 .OR. Vegetation_Type(n) == 13) then
         ice_coverage(n) = min(ice_coverage(n) + land_coverage(n), 1.0_kind_real)
         land_coverage(n) = 0.0_kind_real
      endif

   endif

   if (lwind) then

     !Interpolate lowest level winds to observation location 
     uu5 = ( u_srf (n,1)*w00 + u_srf (n,2)*w10 + u_srf (n,3)*w01 + u_srf (n,4)*w11 ) * dtsfc  + &
           ( u_srfp(n,1)*w00 + u_srfp(n,2)*w10 + u_srfp(n,3)*w01 + u_srfp(n,4)*w11 ) * dtsfcp 
     vv5 = ( v_srf (n,1)*w00 + v_srf (n,2)*w10 + v_srf (n,3)*w01 + v_srf (n,4)*w11 ) * dtsfc  + &
           ( v_srfp(n,1)*w00 + v_srfp(n,2)*w10 + v_srfp(n,3)*w01 + v_srfp(n,4)*w11 ) * dtsfcp

     sfc_speed = f10*sqrt(uu5*uu5+vv5*vv5)
     wind10    = sfc_speed 
     if (uu5*f10 >= 0.0_kind_real .and. vv5*f10 >= 0.0_kind_real) iquadrant = 1 
     if (uu5*f10 >= 0.0_kind_real .and. vv5*f10 <  0.0_kind_real) iquadrant = 2 
     if (uu5*f10 <  0.0_kind_real .and. vv5*f10 >= 0.0_kind_real) iquadrant = 4 
     if (uu5*f10 <  0.0_kind_real .and. vv5*f10 <  0.0_kind_real) iquadrant = 3 
     if (abs(vv5*f10) >= windlimit) then 
         windratio = (uu5*f10) / (vv5*f10) 
     else 
         windratio = 0.0_kind_real 
         if (abs(uu5*f10) > windlimit) then 
             windratio = windscale * uu5*f10 
         endif 
     endif 
     windangle        = atan(abs(windratio))   ! wind azimuth is in radians 
     wind10_direction = quadcof(iquadrant, 1) * pi + windangle * quadcof(iquadrant, 2)   
     wind_speed(n)           = sfc_speed
     wind_direction(n)       = rad2deg*wind10_direction

   else

     wind_speed(n)           = 0.0_kind_real
     wind_direction(n)       = 0.0_kind_real

   endif

   water_temperature(n)     = max(ts(0) + dtskin(0), 270._kind_real)

   !TODO, is nst_gsi ever > 1?
   !if(nst_gsi > 1 .and. water_coverage(1) > 0.0_kind_real) then
      !water_temperature(n)  = max(data_s(itref)+data_s(idtw)-data_s(idtc) + dtskin(0), 271._kind_real)
   !endif

   land_temperature(n)      = ts(1) + dtskin(1)
   ice_temperature(n)       = min(ts(2) + dtskin(2), 280._kind_real)
   snow_temperature(n)      = min(ts(3) + dtskin(3), 280._kind_real)
   soil_moisture_content(n) = sm
   vegetation_fraction(n)   = vfr
   soil_temperature(n)      = stp
   snow_depth(n)            = sn

 enddo

 deallocate(map_to_crtm_ir)
 deallocate(map_to_crtm_mwave)

end subroutine crtm_surface

!----------------------------------------------------------------------------

subroutine kdtree_getweightsindex(comm, geom, nobs, ngrid, lats_ob, lons_ob, nn, interp_w, interp_i)

implicit none

!Arguments
type(fckit_mpi_comm), intent(in)  :: comm               !Communicator
type(fv3jedi_geom),   intent(in)  :: geom               !Model geometry
integer,              intent(in)  :: nobs               !Number of obs on this processor
integer,              intent(in)  :: ngrid              !Number of grid points on this processor
real(kind=kind_real), intent(in)  :: lats_ob(nobs)      !Observation locations, lats
real(kind=kind_real), intent(in)  :: lons_ob(nobs)      !Observation locations, lons
integer,              intent(in)  :: nn                 !Number of neighbours to get back
real(kind=kind_real), intent(out) :: interp_w(nobs,nn)  !Interpolation weights
integer,              intent(out) :: interp_i(nobs,nn)  !Interpolation indices (global unstructured grid)

!Locals
type(mpl_type)    :: mpl
type(kdtree_type) :: kdtree
integer :: i, j, ob, jj, kk, ngrid_loc, ngrid_all
integer, allocatable :: displs(:), rcvcnt(:), nn_index(:)
logical, allocatable :: mask(:)
real(kind=kind_real) :: wprod, dist, bw(nn), bsw
real(kind=kind_real), allocatable :: lat_loc(:), lon_loc(:), lat_all(:), lon_all(:)
real(kind=kind_real), allocatable :: nn_dist(:)


! Gather the model grid to all processors
! ---------------------------------------

! Allocate arrays for unstructured lat/lon
ngrid_loc = ngrid
ngrid_all = 6 * (geom%npx-1) * (geom%npy-1)

allocate(lat_loc(ngrid_loc))
allocate(lon_loc(ngrid_loc))
allocate(lat_all(ngrid_all))
allocate(lon_all(ngrid_all))

! Fill local unstructured lat/lon
jj = 0
do j = geom%jsc,geom%jec
  do i = geom%isc,geom%iec
     jj = jj + 1
     lat_loc(jj) = geom%grid_lat(i,j)
     lon_loc(jj) = geom%grid_lon(i,j) - pi   !BUMP needs -180 to 180
  enddo
enddo

! Gather receive count from each processor
allocate(rcvcnt(comm%size()))
call comm%allgather(jj, rcvcnt)

! Diplacement for each processor
allocate(displs(comm%size()))
displs(1) = 0
do j = 2,comm%size()
   displs(j) = displs(j-1) + rcvcnt(j-1)
enddo

! Gather the lat/lons to all
call comm%allgather(lat_loc,lat_all,jj,rcvcnt,displs)
call comm%allgather(lon_loc,lon_all,jj,rcvcnt,displs)

! Deallocate
deallocate(rcvcnt,displs)
deallocate(lat_loc,lon_loc)


! Create a KDTree for finding nearest neighbours
! ----------------------------------------------

! No mask to true everywhere
allocate(mask(ngrid_all))
mask = .true.

! Initialize mpl (note kdtree not parallel)
call mpl%init()

! Create kdtree
call kdtree%alloc(mpl,ngrid_all,mask)
call kdtree%init(mpl,lon_all,lat_all)


! Loop over observations calling kdtree to generate barycentric interpolation weights
! -----------------------------------------------------------------------------------

allocate(nn_index(nn))
allocate(nn_dist(nn))

do ob = 1,nobs

  nn_index = 0
  nn_dist = 0.0_kind_real

  call kdtree%find_nearest_neighbors(mpl,lons_ob(ob)-pi,lats_ob(ob),nn,nn_index,nn_dist)

  !if (comm%rank()==0) then
  !  print*, '==========================='
  !  print*, rad2deg*lats_ob(ob)
  !  print*, rad2deg*lat_all(nn_index)
  !  print*, ' '
  !  print*, rad2deg*lons_ob(ob)
  !  print*, rad2deg*lon_all(nn_index) + 180.0
  !  print*, '==========================='
  !endif

  !Barycentric weights formula
  bw(:) = 0.0_kind_real
  do jj = 1,nn
    wprod = 1.0_kind_real
    do kk = 1,nn
      if (jj.ne.kk) then
        call sphere_dist(rad2deg*lon_all(nn_index(jj)),rad2deg*lat_all(nn_index(jj)),&
                         rad2deg*lon_all(nn_index(kk)),rad2deg*lat_all(nn_index(kk)),dist)
        wprod = wprod * dist
      endif
    enddo
    bw(jj) = 1.0_kind_real / wprod
  enddo

  !Barycentric weights
  bsw = 0.0_kind_real
  do jj = 1,nn
    bsw = bsw + (bw(jj) / nn_dist(jj))
  enddo

  interp_w(ob,:) = 0.0_kind_real
  do jj = 1,nn
    interp_w(ob,jj) = ( bw(jj) / nn_dist(jj) ) / bsw
    interp_i(ob,jj) = nn_index(jj)
  enddo

enddo

!Deallocate
call kdtree%dealloc
deallocate(mask)
deallocate(nn_index,nn_dist)
deallocate(lon_all,lat_all)

end subroutine kdtree_getweightsindex

!----------------------------------------------------------------------------

subroutine kdtree_getneighbours(comm, geom, nobs, ngrid, nn, interp_i, field_in, field_out)

!Arguments
type(fckit_mpi_comm), intent(in)  :: comm               !Communicator
type(fv3jedi_geom)  , intent(in)  :: geom               !Model geometry
integer             , intent(in)  :: nobs               !Number of obs on this processor
integer             , intent(in)  :: ngrid              !Number of grid points on this processor
integer             , intent(in)  :: nn                 !Number of neighbours
integer             , intent(in)  :: interp_i(nobs,nn)  !Interpolation index
real(kind=kind_real), intent(in)  :: field_in(geom%isc:geom%iec,geom%jsc:geom%jec) !Fields in
real(kind=kind_real), intent(out) :: field_out(nobs,nn)                            !Field nearest neighbours

!Locals
integer :: i, j, jj, ob, ngrid_loc, ngrid_all
integer, allocatable :: displs(:), rcvcnt(:)
real(kind=kind_real), allocatable :: field_loc(:), field_all(:)


! Gather the field to all processors
! ---------------------------------------

! Allocate arrays for unstructured fields
ngrid_loc = ngrid
ngrid_all = 6 * (geom%npx-1) * (geom%npy-1)

allocate(field_loc(ngrid_loc))
allocate(field_all(ngrid_all))

! Fill local unstructured field
jj = 0
do j = geom%jsc,geom%jec
  do i = geom%isc,geom%iec
     jj = jj + 1
     field_loc(jj) = field_in(i,j)
  enddo
enddo

! Gather receive count from each processor
allocate(rcvcnt(comm%size()))
call comm%allgather(jj, rcvcnt)

! Diplacement for each processor
allocate(displs(comm%size()))
displs(1) = 0
do j = 2,comm%size()
   displs(j) = displs(j-1) + rcvcnt(j-1)
enddo

! Gather the lat/lons to all
call comm%allgather(field_loc,field_all,jj,rcvcnt,displs)

! Deallocate
deallocate(rcvcnt,displs)
deallocate(field_loc)


! Set Field out using indices from the kdtree
! -------------------------------------------

do ob = 1,nobs
   field_out(ob,1) = field_all(interp_i(ob,1))
   field_out(ob,2) = field_all(interp_i(ob,2))
   field_out(ob,3) = field_all(interp_i(ob,3))
   field_out(ob,4) = field_all(interp_i(ob,4))
enddo

end subroutine kdtree_getneighbours

!----------------------------------------------------------------------------

subroutine get_lai(lai_type,lai)

  implicit none

  integer             , intent(in ) :: lai_type
  real(kind=kind_real), intent(out) :: lai

  !Dummy code, needs to be figured out
  if (lai_type .ne. 0) then
    lai = 0.0_kind_real
  endif

  end subroutine get_lai

!----------------------------------------------------------------------------

end module surface_vt_mod
