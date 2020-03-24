! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module surface_vt_mod

use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_kinds_mod, only: kind_real
use crtm_module, only: crtm_irlandcoeff_classification
use fv3jedi_constants_mod, only: rad2deg, deg2rad, pi
use fckit_mpi_module, only: fckit_mpi_comm, fckit_mpi_sum
use fckit_geometry_module, only: sphere_distance
use fckit_kdtree_module, only: kdtree,kdtree_create,kdtree_destroy,kdtree_k_nearest_neighbors

use unstructured_interpolation_mod, only: unstrc_interp

implicit none
private

public crtm_surface

contains

!----------------------------------------------------------------------------
! Surface quantities in the form needed by the crtm -------------------------
!----------------------------------------------------------------------------

subroutine crtm_surface( geom, nobs, ngrid, lats_ob, lons_ob, &
                         fld_slmsk, fld_sheleg, fld_tsea, fld_vtype, fld_stype, fld_vfrac, fld_stc, &
                         fld_smc, fld_u_srf, fld_v_srf, fld_f10m, &
                         land_type, vegetation_type, soil_type, water_coverage, land_coverage, ice_coverage, &
                         snow_coverage, lai, water_temperature, land_temperature, ice_temperature, &
                         snow_temperature, soil_moisture_content, vegetation_fraction, soil_temperature, snow_depth, &
                         wind_speed, wind_direction,fld_sss , sea_surface_salinity)

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
real(kind=kind_real), intent(in)  :: fld_stc   (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: fld_smc   (geom%isc:geom%iec,geom%jsc:geom%jec,1)
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
real(kind=kind_real), optional, intent(in)  :: fld_sss   (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), optional, intent(out) :: sea_surface_salinity(nobs)

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
real(kind=kind_real) :: ss00,ss01,ss10,ss11
real(kind=kind_real) :: tsavg,wgtmin,ssavg
real(kind=kind_real) :: vty, sty, vfr, stp, sm, sn, ss
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
real(kind=kind_real), allocatable, dimension(:,:) :: sss

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
real(kind=kind_real), allocatable, dimension(:,:) :: sssp

type(fckit_mpi_comm) :: comm
integer :: i, j, jj
real(kind=kind_real), allocatable, dimension(:) :: slmsk_ob, sheleg_ob, tsea_ob, vtype_ob, stype_ob, &
                                                   vfrac_ob, stc_ob, smc_ob, u_srf_ob, v_srf_ob, f10m_ob, sss_ob
real(kind=kind_real), allocatable :: lat_loc(:), lon_loc(:)
type(unstrc_interp) :: unsinterp

! Communicator from geometry
comm = geom%f_comm

! Number of neighbors
nn = 4

allocate(slmsk(nn,nobs))
allocate(sheleg(nn,nobs))
allocate(tsea(nn,nobs))
allocate(vtype(nn,nobs))
allocate(stype(nn,nobs))
allocate(vfrac(nn,nobs))
allocate(stc(nn,nobs))
allocate(smc(nn,nobs))
allocate(snwdph(nn,nobs))
allocate(u_srf(nn,nobs))
allocate(v_srf(nn,nobs))
allocate(f10m(nn,nobs))
allocate(sss(nn,nobs))

allocate(rslmsk(nn,nobs))
allocate(rvtype(nn,nobs))
allocate(rstype(nn,nobs))

allocate(slmsk_ob (nobs))
allocate(sheleg_ob(nobs))
allocate(tsea_ob  (nobs))
allocate(vtype_ob (nobs))
allocate(stype_ob (nobs))
allocate(vfrac_ob (nobs))
allocate(stc_ob   (nobs))
allocate(smc_ob   (nobs))
allocate(u_srf_ob (nobs))
allocate(v_srf_ob (nobs))
allocate(f10m_ob  (nobs))
allocate(sss_ob   (nobs))

allocate(slmskp(nn,nobs))
allocate(shelegp(nn,nobs))
allocate(tseap(nn,nobs))
allocate(vtypep(nn,nobs))
allocate(stypep(nn,nobs))
allocate(vfracp(nn,nobs))
allocate(stcp(nn,nobs))
allocate(smcp(nn,nobs))
allocate(snwdphp(nn,nobs))
allocate(u_srfp(nn,nobs))
allocate(v_srfp(nn,nobs))
allocate(f10mp(nn,nobs))
allocate(sssp(nn,nobs))

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
 sssp = 0.0_kind_real

 ! Fill local unstructured lat/lon
 allocate(lat_loc(ngrid))
 allocate(lon_loc(ngrid))

 jj = 0
 do j = geom%jsc,geom%jec
   do i = geom%isc,geom%iec
      jj = jj + 1
      lat_loc(jj) = rad2deg*geom%grid_lat(i,j)
      lon_loc(jj) = rad2deg*geom%grid_lon(i,j)
   enddo
 enddo

 call unsinterp%create( comm, nn, 'barycent', ngrid, lat_loc, lon_loc, &
                                              nobs,  lats_ob, lons_ob )

 deallocate(lat_loc,lon_loc)

 ! Get neighbours
 call unsinterp%apply(fld_slmsk (:,:,1), slmsk_ob , rslmsk)
 call unsinterp%apply(fld_sheleg(:,:,1), sheleg_ob, sheleg)
 call unsinterp%apply(fld_tsea  (:,:,1), tsea_ob  , tsea  )
 call unsinterp%apply(fld_vtype (:,:,1), vtype_ob , rvtype)
 call unsinterp%apply(fld_stype (:,:,1), stype_ob , rstype)
 call unsinterp%apply(fld_vfrac (:,:,1), vfrac_ob , vfrac )
 call unsinterp%apply(fld_stc   (:,:,1), stc_ob   , stc   )
 call unsinterp%apply(fld_smc   (:,:,1), smc_ob   , smc   )
 call unsinterp%apply(fld_u_srf (:,:,1), u_srf_ob , u_srf )
 call unsinterp%apply(fld_v_srf (:,:,1), v_srf_ob , v_srf )
 call unsinterp%apply(fld_f10m  (:,:,1), f10m_ob  , f10m  )
 if ( present(fld_sss) ) then
   call unsinterp%apply(fld_sss   (:,:,1), sss_ob  , sss  )
 endif

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

    w00 = unsinterp%interp_w(1,n)
    w10 = unsinterp%interp_w(2,n)
    w01 = unsinterp%interp_w(3,n)
    w11 = unsinterp%interp_w(4,n)

    istyp00 = slmsk(1,n)
    istyp10 = slmsk(2,n)
    istyp01 = slmsk(3,n)
    istyp11 = slmsk(4,n)

    !sno00 = snwdph(1,n)*dtsfc + snwdphp(1,n)*dtsfcp
    !sno01 = snwdph(2,n)*dtsfc + snwdphp(2,n)*dtsfcp
    !sno10 = snwdph(3,n)*dtsfc + snwdphp(3,n)*dtsfcp
    !sno11 = snwdph(4,n)*dtsfc + snwdphp(4,n)*dtsfcp

    sno00 = sheleg(1,n)*dtsfc + shelegp(1,n)*dtsfcp
    sno01 = sheleg(2,n)*dtsfc + shelegp(2,n)*dtsfcp
    sno10 = sheleg(3,n)*dtsfc + shelegp(3,n)*dtsfcp
    sno11 = sheleg(4,n)*dtsfc + shelegp(4,n)*dtsfcp

    sst00 = tsea(1,n)*dtsfc + tsea(1,n)*dtsfcp
    sst01 = tsea(2,n)*dtsfc + tsea(2,n)*dtsfcp
    sst10 = tsea(3,n)*dtsfc + tsea(3,n)*dtsfcp
    sst11 = tsea(4,n)*dtsfc + tsea(4,n)*dtsfcp

    tsavg = sst00*w00 + sst10*w10 + sst01*w01 + sst11*w11

    if ( present(fld_sss) ) then
       ss00 = sss(1,n)*dtsfc + sss(1,n)*dtsfcp
       ss01 = sss(2,n)*dtsfc + sss(2,n)*dtsfcp
       ss10 = sss(3,n)*dtsfc + sss(3,n)*dtsfcp
       ss11 = sss(4,n)*dtsfc + sss(4,n)*dtsfcp

       ssavg = ss00*w00 + ss10*w10 + ss01*w01 + ss11*w11
    endif

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
    ss=0.0_kind_real

    idomsfc=slmsk(1,n)
    wgtmin = w00

    if(istyp00 == 1)then
       vty  = vtype(1,n)
       sty  = stype(1,n)
       wgtavg(1) = wgtavg(1) + w00
       ts(1)=ts(1)+w00*sst00
       vfr  =vfr  +w00*( vfrac(1,n) * dtsfc + &
                        vfracp(1,n) * dtsfcp  )
       stp  =stp  +w00*(   stc(1,n) * dtsfc + &
                          stcp(1,n) * dtsfcp  )
       sm   =sm   +w00*(   smc(1,n) * dtsfc + &
                          smcp(1,n) * dtsfcp  )
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
       if ( present(fld_sss) ) ss   =ss   +w00*ss00
    end if

    if(istyp01 == 1)then
       if(wgtmin < w01 .or. (vty == 0.0_kind_real .and. sty == 0.0_kind_real))then
          vty  = vtype(3,n)
          sty  = stype(3,n)
       end if
       wgtavg(1) = wgtavg(1) + w01
       ts(1)=ts(1)+w01*sst01
       vfr  =vfr  +w01*( vfrac(3,n) * dtsfc + &
                        vfracp(3,n) * dtsfcp  )
       stp  =stp  +w01*(   stc(3,n) * dtsfc + &
                          stcp(3,n) * dtsfcp  )
       sm   =sm   +w01*(   smc(3,n) * dtsfc + &
                          smcp(3,n) * dtsfcp  )
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
       if ( present(fld_sss) ) ss   =ss   +w01*ss01
    end if
    if(wgtmin < w01)then
       idomsfc=slmsk(3,n)
       wgtmin = w01
    end if

    if(istyp10 == 1)then
       if(wgtmin < w10 .or. (vty == 0.0_kind_real .and. sty == 0.0_kind_real))then
          vty  = vtype(2,n)
          sty  = stype(2,n)
       end if
       wgtavg(1) = wgtavg(1) + w10
       ts(1)=ts(1)+w10*sst10
       vfr  =vfr  +w10*(vfrac (2,n) * dtsfc + &
                        vfracp(2,n) * dtsfcp  )
       stp  =stp  +w10*(  stc (2,n) * dtsfc + &
                          stcp(2,n) * dtsfcp  )
       sm   =sm   +w10*(  smc (2,n) * dtsfc + &
                          smcp(2,n) * dtsfcp  )
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
       if ( present(fld_sss) ) ss = ss + w10*ss10
    end if
    if(wgtmin < w10)then
       idomsfc=slmsk(2,n)
       wgtmin = w10
    end if

    if(istyp11 == 1)then
       if(wgtmin < w11 .or. (vty == 0.0_kind_real .and. sty == 0.0_kind_real))then
          vty  = vtype(4,n)
          sty  = stype(4,n)
       endif
       wgtavg(1) = wgtavg(1) + w11
       ts(1)=ts(1)+w11*sst11
       vfr  =vfr  +w11*(vfrac (4,n) * dtsfc + &
                        vfracp(4,n) * dtsfcp  )
       stp  =stp  +w11*(  stc (4,n) * dtsfc + &
                          stcp(4,n) * dtsfcp  )
       sm   =sm   +w11*(  smc (4,n) * dtsfc + &
                          smcp(4,n) * dtsfcp  )
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
       if ( present(fld_sss) ) ss = ss + w11*ss11
    end if

    if(wgtmin < w11)then
       idomsfc=slmsk(4,n)
       wgtmin = w11
    end if

    if(wgtavg(0) > 0.0_kind_real)then
       ts(0) = ts(0)/wgtavg(0)
       if ( present(fld_sss) ) ss    = ss/wgtavg(0)
    else
       ts(0) = tsavg
       if ( present(fld_sss) ) ss    = ssavg
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

    f10 = ( f10m (1,n)*w00 + f10m (2,n)*w10 + f10m (3,n)*w01 + f10m (4,n)*w11 ) * dtsfc + &
          ( f10mp(1,n)*w00 + f10mp(2,n)*w10 + f10mp(3,n)*w01 + f10mp(4,n)*w11 ) * dtsfcp

! Stage 2 - like crtm_interface from GSI
! --------------------------------------

   itype  = vty
   istype = sty

   itype  = min(max(1,itype),nvege_type)
   istype = min(max(1,istype),SOIL_N_TYPES)
   land_type(n) = max(1,map_to_crtm_mwave(itype))
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
     uu5 = ( u_srf (1,n)*w00 + u_srf (2,n)*w10 + u_srf (3,n)*w01 + u_srf (4,n)*w11 ) * dtsfc  + &
           ( u_srfp(1,n)*w00 + u_srfp(2,n)*w10 + u_srfp(3,n)*w01 + u_srfp(4,n)*w11 ) * dtsfcp
     vv5 = ( v_srf (1,n)*w00 + v_srf (2,n)*w10 + v_srf (3,n)*w01 + v_srf (4,n)*w11 ) * dtsfc  + &
           ( v_srfp(1,n)*w00 + v_srfp(2,n)*w10 + v_srfp(3,n)*w01 + v_srfp(4,n)*w11 ) * dtsfcp

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
   if ( present(fld_sss) ) sea_surface_salinity(n)  = ss

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

 call unsinterp%delete()
 deallocate(map_to_crtm_ir)
 deallocate(map_to_crtm_mwave)

end subroutine crtm_surface

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
