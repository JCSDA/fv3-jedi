! (C) Copyright 2018 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!> Variable transforms on wind variables for fv3-jedi 
!> Daniel Holdaway, NASA/JCSDA

module surface_vt_mod

use fv3jedi_geom_mod, only: fv3jedi_geom
use kinds, only: kind_real
use crtm_module, only: crtm_irlandcoeff_classification
use fv3jedi_constants, only: rad2deg, deg2rad, pi

implicit none
public

contains

!----------------------------------------------------------------------------
! Surface quantities in the form needed by the crtm ------------------------- 
!----------------------------------------------------------------------------

subroutine crtm_surface( geom, nobs, nn, weights, slmsk, sheleg, tsea, vtype, & 
                         stype, vfrac, stc, smc, snwdph, u_srf, v_srf, f10m, &
                         land_type, vegetation_type, soil_type, water_coverage, land_coverage, ice_coverage, &
                         snow_coverage, lai, water_temperature, land_temperature, ice_temperature, &
                         snow_temperature, soil_moisture_content, vegetation_fraction, soil_temperature, snow_depth, &
                         wind_speed, wind_direction )

 implicit none

 !Arguments
 type(fv3jedi_geom)  , intent(in ) :: geom !Geometry for the model
 integer, intent(in)               :: nobs, nn

 real(kind=kind_real), intent(inout), dimension(nobs,nn) :: weights    !Weights applied to neighbours
 integer             , intent(inout), dimension(nobs,nn) :: slmsk
 real(kind=kind_real), intent(inout), dimension(nobs,nn) :: sheleg
 real(kind=kind_real), intent(inout), dimension(nobs,nn) :: tsea
 integer             , intent(inout), dimension(nobs,nn) :: vtype
 integer             , intent(inout), dimension(nobs,nn) :: stype
 real(kind=kind_real), intent(inout), dimension(nobs,nn) :: vfrac
 real(kind=kind_real), intent(inout), dimension(nobs,nn) :: snwdph
 real(kind=kind_real), intent(inout), dimension(nobs,nn) :: stc
 real(kind=kind_real), intent(inout), dimension(nobs,nn) :: smc
 real(kind=kind_real), intent(inout), dimension(nobs,nn) :: u_srf
 real(kind=kind_real), intent(inout), dimension(nobs,nn) :: v_srf
 real(kind=kind_real), intent(inout), dimension(nobs,nn) :: f10m

 !Incase of need for two time levels
 integer             , dimension(nobs,nn) :: slmskp
 real(kind=kind_real), dimension(nobs,nn) :: shelegp
 real(kind=kind_real), dimension(nobs,nn) :: tseap
 integer             , dimension(nobs,nn) :: vtypep
 integer             , dimension(nobs,nn) :: stypep
 real(kind=kind_real), dimension(nobs,nn) :: vfracp
 real(kind=kind_real), dimension(nobs,nn) :: snwdphp
 real(kind=kind_real), dimension(nobs,nn) :: stcp
 real(kind=kind_real), dimension(nobs,nn) :: smcp
 real(kind=kind_real), dimension(nobs,nn) :: u_srfp
 real(kind=kind_real), dimension(nobs,nn) :: v_srfp
 real(kind=kind_real), dimension(nobs,nn) :: f10mp

 integer             , intent(out), dimension(nobs) :: vegetation_type
 integer             , intent(out), dimension(nobs) :: land_type
 integer             , intent(out), dimension(nobs) :: soil_type
 real(kind=kind_real), intent(out), dimension(nobs) :: water_coverage
 real(kind=kind_real), intent(out), dimension(nobs) :: land_coverage
 real(kind=kind_real), intent(out), dimension(nobs) :: ice_coverage
 real(kind=kind_real), intent(out), dimension(nobs) :: snow_coverage
 real(kind=kind_real), intent(out), dimension(nobs) :: lai
 real(kind=kind_real), intent(out), dimension(nobs) :: water_temperature
 real(kind=kind_real), intent(out), dimension(nobs) :: land_temperature
 real(kind=kind_real), intent(out), dimension(nobs) :: ice_temperature
 real(kind=kind_real), intent(out), dimension(nobs) :: snow_temperature
 real(kind=kind_real), intent(out), dimension(nobs) :: soil_moisture_content
 real(kind=kind_real), intent(out), dimension(nobs) :: vegetation_fraction
 real(kind=kind_real), intent(out), dimension(nobs) :: soil_temperature
 real(kind=kind_real), intent(out), dimension(nobs) :: snow_depth
 real(kind=kind_real), intent(out), dimension(nobs) :: wind_speed
 real(kind=kind_real), intent(out), dimension(nobs) :: wind_direction

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


 !Second time level option, zero for now
 slmskp = 0.0_kind_real
 shelegp = 0.0_kind_real
 tseap = 0.0_kind_real
 vtypep = 0.0_kind_real
 stypep = 0.0_kind_real
 vfracp = 0.0_kind_real
 snwdphp = 0.0_kind_real
 stcp = 0.0_kind_real
 smcp = 0.0_kind_real

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

    w00 = weights(n,1)
    w10 = weights(n,2)
    w01 = weights(n,3)
    w11 = weights(n,4)
 
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

subroutine crtm_surface_neighbours( geom, nobs, ngrid, lats_ob, lons_ob, nn, weights, &
                                    slmsk, sheleg, tsea, vtype, stype, vfrac, stc, smc, snwdph, &
                                    u_srf, v_srf, f10m, &
                                    slmsko, shelego, tseao, vtypeo, stypeo, vfraco, stco, smco, snwdpho, &
                                    u_srfo, v_srfo, f10mo )

use mpp_mod, only: mpp_npes, mpp_pe
use mpi

!Return the neighbouring grid points and weights 

implicit none

type(fv3jedi_geom)  , intent(in)    :: geom               !Model geometry
integer             , intent(in)    :: nobs               !Number of obs on this processor
integer             , intent(in)    :: ngrid              !Number of obs on this processor
real(kind=kind_real), intent(in)    :: lats_ob(nobs)      !Observation locations, lats
real(kind=kind_real), intent(in)    :: lons_ob(nobs)      !Observation locations, lons
integer             , intent(in)    :: nn                 !Number of neighbours

integer             , intent(in), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: slmsk
real(kind=kind_real), intent(in), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: sheleg
real(kind=kind_real), intent(in), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: tsea
integer             , intent(in), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: vtype
integer             , intent(in), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: stype
real(kind=kind_real), intent(in), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: vfrac
real(kind=kind_real), intent(in), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: snwdph
real(kind=kind_real), intent(in), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,4) :: stc
real(kind=kind_real), intent(in), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed,4) :: smc
real(kind=kind_real), intent(in), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: u_srf
real(kind=kind_real), intent(in), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: v_srf
real(kind=kind_real), intent(in), dimension(geom%bd%isd:geom%bd%ied,geom%bd%jsd:geom%bd%jed)   :: f10m

!Outputs on stencil used
real(kind=kind_real), intent(inout), dimension(nobs,nn) :: weights    !Weights applied to neighbours
integer             , intent(inout), dimension(nobs,nn) :: slmsko
real(kind=kind_real), intent(inout), dimension(nobs,nn) :: shelego
real(kind=kind_real), intent(inout), dimension(nobs,nn) :: tseao
integer             , intent(inout), dimension(nobs,nn) :: vtypeo
integer             , intent(inout), dimension(nobs,nn) :: stypeo
real(kind=kind_real), intent(inout), dimension(nobs,nn) :: vfraco
real(kind=kind_real), intent(inout), dimension(nobs,nn) :: snwdpho
real(kind=kind_real), intent(inout), dimension(nobs,nn) :: stco
real(kind=kind_real), intent(inout), dimension(nobs,nn) :: smco
real(kind=kind_real), intent(inout), dimension(nobs,nn) :: u_srfo
real(kind=kind_real), intent(inout), dimension(nobs,nn) :: v_srfo
real(kind=kind_real), intent(inout), dimension(nobs,nn) :: f10mo

real(kind=kind_real), allocatable, dimension(:,:) :: grid_loc, grid_glo
real(kind=kind_real), allocatable, dimension(:)   :: grid_tmp
real(kind=kind_real), allocatable, dimension(:,:,:,:) :: grid_glo_str
real(kind=kind_real), allocatable, dimension(:,:,:) :: grid_obs

integer :: isc,iec,jsc,jec, nvars, ngrid_glo, lowerb, upperb
integer :: i,j,jj,n,t
integer :: tind, iind, jind

integer :: npes, peid, ierr, stquad
integer, allocatable, dimension(:) :: ngridv, displs, rcvcnt

real(kind=kind_real) :: llo(2), ll1(2), ll2(2), ll3(2), ll4(2), dist(4)
real(kind=kind_real) :: dist_new, lats_gd, lons_gd

integer :: ind_cor, ind_til, ind_ict, ind_jct, ind_lat, ind_lon, ind_wei
integer :: tmin(4), imin(4), jmin(4)


integer :: i1, i2, j1, j2

! Processor info
! --------------
npes = mpp_npes()
peid = mpp_pe()


! Convenience
! -----------
isc = geom%bd%isc
iec = geom%bd%iec
jsc = geom%bd%jsc
jec = geom%bd%jec


! Create array holding all information that needs to be passed to ob location
! ---------------------------------------------------------------------------

nvars = 7 !Number of things to be passed

ind_cor = 1
ind_til = 2
ind_ict = 3
ind_jct = 4
ind_lat = 5
ind_lon = 6
ind_wei = 7 !Always last as not allocated until ob procs

allocate(grid_loc(ngrid,nvars-1))

jj = 0
do j = jsc,jec
  do i = isc,iec

     jj = jj + 1

     !Bookkeeping of corners as the number of neighbours is different
     grid_loc(jj,ind_cor) = 0.0_kind_real
     if (i == 1          .and. j == 1         ) grid_loc(jj,ind_cor) = 1.0_kind_real !SW Corner
     if (i == geom%npx-1 .and. j == 1         ) grid_loc(jj,ind_cor) = 2.0_kind_real !SE Corner
     if (i == geom%npx-1 .and. j == geom%npy-1) grid_loc(jj,ind_cor) = 3.0_kind_real !NE Corner
     if (i == 1          .and. j == geom%npy-1) grid_loc(jj,ind_cor) = 4.0_kind_real !NW Corner

     grid_loc(jj,ind_til) = real(geom%ntile,kind_real)
     grid_loc(jj,ind_ict) = real(i,kind_real)
     grid_loc(jj,ind_jct) = real(j,kind_real)
     grid_loc(jj,ind_lat) = deg2rad*geom%grid_lat(i,j)
     grid_loc(jj,ind_lon) = deg2rad*geom%grid_lon(i,j)

  enddo
enddo

! Gather the model grid to all processors
! ---------------------------------------

allocate(ngridv(0:npes-1))
allocate(displs(0:npes-1))
allocate(rcvcnt(0:npes-1))

ngridv = 0
ngridv(peid) = ngrid

do n = 0,npes-1
   displs(n) = n
   rcvcnt(n) = 1
enddo

call mpi_allgatherv(ngridv(peid), 1, mpi_int, ngridv, rcvcnt, displs, mpi_int, mpi_comm_world, ierr)
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

allocate(grid_glo(ngrid_glo,nvars-1))
allocate(grid_tmp(ngrid_glo))

do n = 1,nvars-1
  grid_tmp(lowerb:upperb) = grid_loc(:,n)
  call mpi_allgatherv(grid_tmp(lowerb), rcvcnt(peid), mpi_real8, grid_tmp, rcvcnt, displs, mpi_real8, mpi_comm_world, ierr)
  grid_glo(:,n) = grid_tmp
enddo

deallocate(ngridv,displs,rcvcnt)
deallocate(grid_tmp,grid_loc)

!Everyone has entire grid_loc


! Reorganise to a structured global grid
! --------------------------------------
allocate(grid_glo_str(nvars-1,geom%ntiles,isc:iec,jsc:jec))
grid_glo_str = 0.0_kind_real

do n = 1,nvars-1
  do j = 1,ngrid_glo

     tind = nint(grid_glo(j,ind_til))
     iind = nint(grid_glo(j,ind_ict))
     jind = nint(grid_glo(j,ind_jct))

     grid_glo_str(n,tind,iind,jind) = grid_glo(j,n)

  enddo
enddo

deallocate(grid_glo)


! Loop over observations and compute weights and create array with fields at the neighbours
! -----------------------------------------------------------------------------------------

allocate(grid_obs(nobs,4,nvars))

do n = 1,nobs

  llo(1) = lons_ob(n)
  llo(2) = lats_ob(n)

  dist = 10.0e10
  do j = jsc,jec
    do i = isc,iec
      do t = 1,geom%ntiles
     
        ll1(1) = grid_glo_str(ind_lon,t,i,j)
        ll1(2) = grid_glo_str(ind_lat,t,i,j)
        call great_circle_dist(llo,ll1,dist_new)

        if (dist_new < dist(1)) then
          dist(4) = dist(3)
          dist(3) = dist(2)
          dist(2) = dist(1)
          dist(1) = dist_new
          tmin(4) = tmin(3)
          tmin(3) = tmin(2)
          tmin(2) = tmin(1)
          tmin(1) = t
          imin(4) = imin(3)
          imin(3) = imin(2)
          imin(2) = imin(1)
          imin(1) = i
          jmin(4) = jmin(3)
          jmin(3) = jmin(2)
          jmin(2) = jmin(1)
          jmin(1) = j
        endif

      enddo
    enddo 
  enddo

  ll1(1) = grid_glo_str(ind_lon,tmin(1),imin(1),jmin(1))
  ll1(2) = grid_glo_str(ind_lat,tmin(1),imin(1),jmin(1))
  ll2(1) = grid_glo_str(ind_lon,tmin(2),imin(2),jmin(2))
  ll2(2) = grid_glo_str(ind_lat,tmin(2),imin(2),jmin(2))
  ll3(1) = grid_glo_str(ind_lon,tmin(3),imin(3),jmin(3))
  ll3(2) = grid_glo_str(ind_lat,tmin(3),imin(3),jmin(3))
  ll4(1) = grid_glo_str(ind_lon,tmin(4),imin(4),jmin(4))
  ll4(2) = grid_glo_str(ind_lat,tmin(4),imin(4),jmin(4))

  call bilin_interp_irregular(ll1,ll2,ll3,ll4,llo,grid_obs(n,:,ind_wei))

 
  !print*, 'corner', grid_glo_str(ind_cor,tmin,imin,jmin)

  !print*, 'lats', grid_glo_str(ind_lat,tmin,i1,j1), grid_glo_str(ind_lat,tmin,i1,j2), grid_glo_str(ind_lat,tmin,i2,j1), grid_glo_str(ind_lat,tmin,i2,j2)
  !print*, 'lons', grid_glo_str(ind_lat,tmin,i1,j1), grid_glo_str(ind_lat,tmin,i1,j2), grid_glo_str(ind_lat,tmin,i2,j1), grid_glo_str(ind_lat,tmin,i2,j2)

call abor1_ftn("crtm_surface_neighbours done")


enddo

deallocate(grid_glo_str)

call abor1_ftn("crtm_surface_neighbours done")

end subroutine crtm_surface_neighbours

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

  subroutine bilin_interp_irregular(ll1,ll2,ll3,ll4,llo,w)

   implicit none
   real(kind=kind_real), intent(in ) :: ll1(2), ll2(2), ll3(2), ll4(2), llo(2)
   real(kind=kind_real), intent(out) :: w(4)

   integer :: indin(4), indout(4),  tind

   real(kind=kind_real) :: lats(4), latssort(2,4), tmp1, tmp2
   real(kind=kind_real) :: lons(4), lonssort(2,4)
   real(kind=kind_real) :: lat1, lat2, lat3, lat4
   real(kind=kind_real) :: lon1, lon2, lon3, lon4

   real(kind=kind_real) :: dll1(2), dll2(2)

   real(kind=kind_real) :: dlat12, dlon12
   real(kind=kind_real) :: dlat23, dlon23
   real(kind=kind_real) :: dlat34, dlon34
   real(kind=kind_real) :: dlat41, dlon41

   real(kind=kind_real) :: s, t, t1, t2

   integer :: j,k

   !Index 1 of ll is lon
   !Index 2 of ll is lat

   indin(1) = 1
   indin(2) = 2
   indin(3) = 3
   indin(4) = 4

   lats(1) = ll1(1)
   lats(2) = ll2(1)
   lats(3) = ll3(1)
   lats(4) = ll4(1)

   lons(1) = ll1(2)
   lons(2) = ll2(2)
   lons(3) = ll3(2)
   lons(4) = ll4(2)

   !Sort to structure the points
   ! 1 = SW
   ! 2 = SE
   ! 3 = NE
   ! 4 = NW
   latssort(1,:) = lats
   latssort(2,:) = lons
   indout = indin
   do j = 1,3
     do k = j,4
       if (latssort(1,j) > latssort(1,k)) then
          tmp1 = latssort(1,k)
          tmp2 = latssort(2,k)
          tind = indout(k)
          latssort(1,k) = latssort(1,j)
          latssort(1,j) = tmp1
          latssort(2,k) = latssort(2,j)
          latssort(2,j) = tmp2
          indout(k) = indout(j)
          indout(j) = tind
       endif
     enddo
   enddo

   if (latssort(2,1) < latssort(2,2)) then
     lat1 = latssort(1,1)
     lon1 = latssort(2,1)
     lat2 = latssort(1,2)
     lon2 = latssort(2,2)
   else
     lat2 = latssort(1,1)
     lon2 = latssort(2,1)
     lat1 = latssort(1,2)
     lon1 = latssort(2,2)
   endif
   if (latssort(2,3) < latssort(2,4)) then
     lat4 = latssort(1,3)
     lon4 = latssort(2,3)
     lat3 = latssort(1,4)
     lon3 = latssort(2,4)
   else
     lat3 = latssort(1,3)
     lon3 = latssort(2,3)
     lat4 = latssort(1,4)
     lon4 = latssort(2,4)
   endif

   print*, '1', lat1, lon1
   print*, '2', lat2, lon2
   print*, '3', lat3, lon3
   print*, '4', lat4, lon4

   dll1(1) = lat1
   dll1(2) = lon1
   dll2(1) = lat2
   dll2(2) = lon1
   call great_circle_dist(dll1,dll2,dlat12)

   dll1(1) = lat1
   dll1(2) = lon1
   dll2(1) = lat1
   dll2(2) = lon2
   call great_circle_dist(dll1,dll2,dlon12)

   dll1(1) = lat2
   dll1(2) = lon2
   dll2(1) = lat3
   dll2(2) = lon2
   call great_circle_dist(dll1,dll2,dlat23)

   dll1(1) = lat2
   dll1(2) = lon2
   dll2(1) = lat2
   dll2(2) = lon3
   call great_circle_dist(dll1,dll2,dlon23)

   dll1(1) = lat3
   dll1(2) = lon3
   dll2(1) = lat4
   dll2(2) = lon3
   call great_circle_dist(dll1,dll2,dlat34)

   dll1(1) = lat3
   dll1(2) = lon3
   dll2(1) = lat3
   dll2(2) = lon4
   call great_circle_dist(dll1,dll2,dlon34)

   dll1(1) = lat4
   dll1(2) = lon4
   dll2(1) = lat1
   dll2(2) = lon4
   call great_circle_dist(dll1,dll2,dlat41)

   dll1(1) = lat4
   dll1(2) = lon4
   dll2(1) = lat4
   dll2(2) = lon1
   call great_circle_dist(dll1,dll2,dlon41)

   A = dlon41 * dlat34 - dlat31 * dlon42
   B = llo(1) * (dlon23 - dlon41) - llo(2) * (dlat23 - dlat41) + dlat41*lat2 - dlat41*lon2 + lat4*dlat23 - lat4*dlon23
   C = llo(1)*dlon34 - llo(2)*dlat34 + lon4*lat3 - lon3*lat4

   t1 = -B - sqrt(B**2 - 4*A*C) / (2*C)
   t2 = -B + sqrt(B**2 - 4*A*C) / (2*C)

   if (0.0 <= t1 .and. t1 <= 1.0) then
      t = t1
   elseif (0.0 <= t2 .and. t2 <= 1.0) then
      t = t2
   else
      call abor1_ftn("bilin_interp_irregular failure")      
   endif

   !check
   if (0.0 <= t1 .and. t1 <= 1.0 .and. 0.0 <= t2 .and. t2 <= 1.0) then
      call abor1_ftn("bilin_interp_irregular no unique solution")
   endif

   s = (llo(1) - lat4 - dlat41*t) / (lat3 + dlat23*t - lat4 - dlat41*t)

   w(1) = (1 - s) * t
   w(2) = s * t
   w(3) = s * (1 - t)
   w(4) = (1 - s) * (1 - t)

  end subroutine bilin_interp_irregular

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

  subroutine great_circle_dist(ll1,ll2,d)

   use fv3jedi_constants, only: radius

   implicit none

   !Arguments
   real(kind=kind_real), intent(in ) :: ll1(2), ll2(2)
   real(kind=kind_real), intent(out) :: d

   !Locals
   real(kind=kind_real) :: dlat, dlon

   dlat = ll2(1) - ll1(1)
   dlon = ll2(2) - ll1(2)

   d = radius * 2 * asin( sqrt( sin(dlat/2.)**2 + cos(ll1(1)) * cos(ll2(1)) * sin(dlon/2.)**2 ) )

  end subroutine great_circle_dist

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------

end module surface_vt_mod
