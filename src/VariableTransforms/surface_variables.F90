! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Variable transforms/interpolation for surface variables in fv3-jedi 
!> Daniel Holdaway, NASA/JCSDA

module surface_vt_mod

use fv3jedi_geom_mod, only: fv3jedi_geom
use kinds, only: kind_real
use crtm_module, only: crtm_irlandcoeff_classification
use fv3jedi_constants, only: rad2deg, deg2rad, pi

implicit none
public crtm_surface

interface crtm_surface_kdtree_getfieldneighbours
 module procedure crtm_surface_kdtree_getfieldneighbours_int
 module procedure crtm_surface_kdtree_getfieldneighbours_real
end interface

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
 integer             , intent(in)  :: fld_slmsk (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)
 real(kind=kind_real), intent(in)  :: fld_sheleg(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)
 real(kind=kind_real), intent(in)  :: fld_tsea  (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)
 integer             , intent(in)  :: fld_vtype (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)
 integer             , intent(in)  :: fld_stype (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)
 real(kind=kind_real), intent(in)  :: fld_vfrac (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)
 real(kind=kind_real), intent(in)  :: fld_stc   (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,4)
 real(kind=kind_real), intent(in)  :: fld_smc   (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,4)
 real(kind=kind_real), intent(in)  :: fld_snwdph(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)
 real(kind=kind_real), intent(in)  :: fld_u_srf (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)
 real(kind=kind_real), intent(in)  :: fld_v_srf (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)
 real(kind=kind_real), intent(in)  :: fld_f10m  (geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)
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

!Number of weights
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

 !Get interpolation weight and index in global grid
 call crtm_surface_kdtree_setup( geom, nobs, ngrid, lats_ob, lons_ob, nn, interp_w, interp_i ) 

 !Get field at nn neighbours
 call crtm_surface_kdtree_getfieldneighbours( geom, nobs, ngrid, nn, interp_i, fld_slmsk     , slmsk  )
 call crtm_surface_kdtree_getfieldneighbours( geom, nobs, ngrid, nn, interp_i, fld_sheleg    , sheleg )
 call crtm_surface_kdtree_getfieldneighbours( geom, nobs, ngrid, nn, interp_i, fld_tsea      , tsea   )
 call crtm_surface_kdtree_getfieldneighbours( geom, nobs, ngrid, nn, interp_i, fld_vtype     , vtype  )
 call crtm_surface_kdtree_getfieldneighbours( geom, nobs, ngrid, nn, interp_i, fld_stype     , stype  )
 call crtm_surface_kdtree_getfieldneighbours( geom, nobs, ngrid, nn, interp_i, fld_vfrac     , vfrac  )
 call crtm_surface_kdtree_getfieldneighbours( geom, nobs, ngrid, nn, interp_i, fld_stc(:,:,1), stc    )
 call crtm_surface_kdtree_getfieldneighbours( geom, nobs, ngrid, nn, interp_i, fld_smc(:,:,1), smc    )
 call crtm_surface_kdtree_getfieldneighbours( geom, nobs, ngrid, nn, interp_i, fld_snwdph    , snwdph )
 call crtm_surface_kdtree_getfieldneighbours( geom, nobs, ngrid, nn, interp_i, fld_u_srf     , u_srf  )
 call crtm_surface_kdtree_getfieldneighbours( geom, nobs, ngrid, nn, interp_i, fld_v_srf     , v_srf  )
 call crtm_surface_kdtree_getfieldneighbours( geom, nobs, ngrid, nn, interp_i, fld_f10m      , f10m   )

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
!----------------------------------------------------------------------------

subroutine crtm_surface_kdtree_setup( geom, nobs, ngrid, lats_ob, lons_ob, nn, interp_w, interp_i )

use mpp_mod, only: mpp_npes, mpp_pe
use mpi
use type_kdtree, only: kdtree_type

implicit none

!Arguments
type(fv3jedi_geom)  , intent(in)    :: geom               !Model geometry
integer             , intent(in)    :: nobs               !Number of obs on this processor
integer             , intent(in)    :: ngrid              !Number of grid points on this processor
real(kind=kind_real), intent(in)    :: lats_ob(nobs)      !Observation locations, lats
real(kind=kind_real), intent(in)    :: lons_ob(nobs)      !Observation locations, lons
integer             , intent(in)    :: nn                 !Number of neighbours to get back
real(kind=kind_real), intent(out)   :: interp_w(nobs,nn)  !Interpolation weights
integer             , intent(out)   :: interp_i(nobs,nn)  !Interpolation indices (global unstructured grid)

!Locals
integer :: npes, peid, ierr
integer :: i, j, k, l, n, jj, lowerb, upperb, ngrid_glo
real(kind=kind_real), allocatable :: grid_lat_loc(:), grid_lon_loc(:), grid_lat_glo(:), grid_lon_glo(:)
integer, allocatable :: ngridv(:), displs(:), rcvcnt(:)
logical, allocatable :: mask(:)
integer, allocatable :: nn_index(:,:)
real(kind=kind_real), allocatable :: nn_dist(:,:)
type(kdtree_type) :: kdtree
real(kind=kind_real) :: dist, tmplat(1), tmplon(1)

! Gather the model grid to all processors
! ---------------------------------------

npes = mpp_npes()
peid = mpp_pe()

!Unstructured local grid
allocate(grid_lat_loc(ngrid))
allocate(grid_lon_loc(ngrid))
jj = 0
do j = geom%bd%jsc,geom%bd%jec
  do i = geom%bd%isc,geom%bd%iec
     jj = jj + 1
     grid_lat_loc(jj) = geom%grid_lat(i,j)
     grid_lon_loc(jj) = geom%grid_lon(i,j)
  enddo
enddo

allocate(ngridv(0:npes-1))
allocate(displs(0:npes-1))
allocate(rcvcnt(0:npes-1))

ngridv = 0
ngridv(peid) = ngrid

do n = 0,npes-1
   displs(n) = n
   rcvcnt(n) = 1
enddo

call mpi_allgatherv(ngrid, 1, mpi_int, ngridv, rcvcnt, displs, mpi_int, mpi_comm_world, ierr)
ngrid_glo = sum(ngridv)

do n = 0,npes-1
   rcvcnt(n) = ngridv(n)
enddo
displs(0) = 0
do n = 1,npes-1
   displs(n) = displs(n-1) + rcvcnt(n-1)
enddo

lowerb = displs(peid)+1
upperb = displs(peid)+rcvcnt(peid)

allocate(grid_lat_glo(ngrid_glo))
allocate(grid_lon_glo(ngrid_glo))

grid_lat_glo(lowerb:upperb) = grid_lat_loc
grid_lon_glo(lowerb:upperb) = grid_lon_loc

call mpi_allgatherv(grid_lat_loc, rcvcnt(peid), mpi_real8, grid_lat_glo, rcvcnt, displs, mpi_real8, mpi_comm_world, ierr)
call mpi_allgatherv(grid_lon_loc, rcvcnt(peid), mpi_real8, grid_lon_glo, rcvcnt, displs, mpi_real8, mpi_comm_world, ierr)

!Create cover tree for gloval grid
allocate(mask(ngrid_glo))
mask = .true.

call kdtree%create(ngrid_glo,deg2rad*grid_lon_glo,deg2rad*grid_lat_glo,mask)


allocate(nn_index(nobs,nn))
allocate(nn_dist(nobs,nn))

do k = 1,nobs

  tmplon(1) = deg2rad*lons_ob(k)
  tmplat(1) = deg2rad*lats_ob(k)

  call kdtree%find_nearest_neighbors(tmplon(1),tmplat(1),nn,nn_index(k,:),nn_dist(k,:))

  dist=sum(nn_dist(k,:))

  do l = 1, nn
     interp_w(k,l) = (dist-nn_dist(k,l))
     interp_i(k,l) = nn_index(k,l)
  end do
  
  interp_w(k,:) = interp_w(k,:) / sum(interp_w(k,:))
  
enddo

!Deallocate
call kdtree%dealloc
deallocate(mask)
deallocate(nn_index,nn_dist)
deallocate(grid_lon_glo,grid_lat_glo)
deallocate(grid_lat_loc,grid_lon_loc)
deallocate(ngridv,displs,rcvcnt)

end subroutine crtm_surface_kdtree_setup

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine crtm_surface_kdtree_getfieldneighbours_int( geom, nobs, ngrid, nn, interp_i, field_in, field_out )

use mpp_mod, only: mpp_npes, mpp_pe
use mpi

implicit none

!Arguments
type(fv3jedi_geom)  , intent(in)  :: geom               !Model geometry
integer             , intent(in)  :: nobs               !Number of obs on this processor
integer             , intent(in)  :: ngrid              !Number of grid points on this processor
integer             , intent(in)  :: nn                 !Number of neighbours
integer             , intent(in)  :: interp_i(nobs,nn)  !Interpolation index

integer             , intent(in)  :: field_in(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed) !Fields in
integer             , intent(out) :: field_out(nobs,nn)                                        !Field nearest neighbours

!Locals
integer :: npes, peid, ierr
integer :: i, j, k, l, n, jj, lowerb, upperb, ngrid_glo
integer, allocatable :: field_loc(:), field_glo(:)
integer, allocatable :: ngridv(:), displs(:), rcvcnt(:)


! Gather the model grid to all processors
! ---------------------------------------

npes = mpp_npes()
peid = mpp_pe()

!Unstructured local grid
allocate(field_loc(ngrid))
jj = 0
do j = geom%bd%jsc,geom%bd%jec
  do i = geom%bd%isc,geom%bd%iec
     jj = jj + 1
     field_loc(jj) = field_in(i,j)
  enddo
enddo

allocate(ngridv(0:npes-1))
allocate(displs(0:npes-1))
allocate(rcvcnt(0:npes-1))

ngridv = 0
ngridv(peid) = ngrid

do n = 0,npes-1
   displs(n) = n
   rcvcnt(n) = 1
enddo

call mpi_allgatherv(ngrid, 1, mpi_int, ngridv, rcvcnt, displs, mpi_int, mpi_comm_world, ierr)
ngrid_glo = sum(ngridv)

do n = 0,npes-1
   rcvcnt(n) = ngridv(n)
enddo
displs(0) = 0
do n = 1,npes-1
   displs(n) = displs(n-1) + rcvcnt(n-1)
enddo

lowerb = displs(peid)+1
upperb = displs(peid)+rcvcnt(peid)

allocate(field_glo(ngrid_glo))

field_glo(lowerb:upperb) = field_loc

call mpi_allgatherv(field_loc, rcvcnt(peid), mpi_int, field_glo, rcvcnt, displs, mpi_int, mpi_comm_world, ierr)


do k = 1,nobs
   field_out(k,1) = field_glo(interp_i(k,1))
   field_out(k,2) = field_glo(interp_i(k,2))
   field_out(k,3) = field_glo(interp_i(k,3))
   field_out(k,4) = field_glo(interp_i(k,4))
enddo


deallocate(field_glo)
deallocate(ngridv,displs,rcvcnt)
deallocate(field_loc)

end subroutine crtm_surface_kdtree_getfieldneighbours_int

subroutine crtm_surface_kdtree_getfieldneighbours_real( geom, nobs, ngrid, nn, interp_i, field_in, field_out )

use mpp_mod, only: mpp_npes, mpp_pe
use mpi

!Arguments
type(fv3jedi_geom)  , intent(in)  :: geom               !Model geometry
integer             , intent(in)  :: nobs               !Number of obs on this processor
integer             , intent(in)  :: ngrid              !Number of grid points on this processor
integer             , intent(in)  :: nn                 !Number of neighbours
integer             , intent(in)  :: interp_i(nobs,nn)  !Interpolation index

real(kind=kind_real), intent(in)  :: field_in(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed) !Fields in
real(kind=kind_real), intent(out) :: field_out(nobs,nn)                                        !Field nearest neighbours

!Locals
integer :: npes, peid, ierr
integer :: i, j, k, l, n, jj, lowerb, upperb, ngrid_glo
real(kind=kind_real), allocatable :: field_loc(:), field_glo(:)
integer, allocatable :: ngridv(:), displs(:), rcvcnt(:)

! Gather the model grid to all processors
! ---------------------------------------

npes = mpp_npes()
peid = mpp_pe()

!Unstructured local grid
allocate(field_loc(ngrid))
jj = 0
do j = geom%bd%jsc,geom%bd%jec
  do i = geom%bd%isc,geom%bd%iec
     jj = jj + 1
     field_loc(jj) = field_in(i,j)
  enddo
enddo

allocate(ngridv(0:npes-1))
allocate(displs(0:npes-1))
allocate(rcvcnt(0:npes-1))

ngridv = 0
ngridv(peid) = ngrid

do n = 0,npes-1
   displs(n) = n
   rcvcnt(n) = 1
enddo

call mpi_allgatherv(ngrid, 1, mpi_int, ngridv, rcvcnt, displs, mpi_int, mpi_comm_world, ierr)
ngrid_glo = sum(ngridv)

do n = 0,npes-1
   rcvcnt(n) = ngridv(n)
enddo
displs(0) = 0
do n = 1,npes-1
   displs(n) = displs(n-1) + rcvcnt(n-1)
enddo

lowerb = displs(peid)+1
upperb = displs(peid)+rcvcnt(peid)

allocate(field_glo(ngrid_glo))

field_glo(lowerb:upperb) = field_loc

call mpi_allgatherv(field_loc, rcvcnt(peid), mpi_real8, field_glo, rcvcnt, displs, mpi_real8, mpi_comm_world, ierr)

do k = 1,nobs
   field_out(k,1) = field_glo(interp_i(k,1))
   field_out(k,2) = field_glo(interp_i(k,2))
   field_out(k,3) = field_glo(interp_i(k,3))
   field_out(k,4) = field_glo(interp_i(k,4))
enddo

deallocate(field_glo)
deallocate(ngridv,displs,rcvcnt)
deallocate(field_loc)

end subroutine crtm_surface_kdtree_getfieldneighbours_real

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine get_lai(lai_type,lai)

  implicit none

  integer             , intent(in ) :: lai_type
  real(kind=kind_real), intent(out) :: lai

  !Needs to be figured out

  end subroutine get_lai

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

end module surface_vt_mod
