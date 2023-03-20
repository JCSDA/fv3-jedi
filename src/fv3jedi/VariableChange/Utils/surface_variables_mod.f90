! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module surface_vt_mod

use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_kinds_mod, only: kind_real
use crtm_module, only: crtm_irlandcoeff_classification
use fv3jedi_constants_mod, only: constant

use slow_unstructured_interpolation_mod, only: unstrc_interp

implicit none
private

public crtm_surface

contains

!----------------------------------------------------------------------------
! Surface quantities in the form needed by the crtm -------------------------
!----------------------------------------------------------------------------

subroutine crtm_surface( geom, day_of_year, &
                         field_slmsk, field_sheleg, field_tsea, field_vtype, field_stype, &
                         field_vfrac, field_stc, field_smc, field_u_srf, field_v_srf, field_f10m, &
                         field_sss, land_type_npoess, land_type_igbp, &
                         vegetation_type, soil_type, water_coverage, land_coverage, &
                         ice_coverage, snow_coverage, lai, water_temperature, land_temperature, &
                         ice_temperature, snow_temperature, soil_moisture_content, &
                         vegetation_fraction, soil_temperature, snow_depth, wind_speed, &
                         wind_direction, sea_surface_salinity)

implicit none

!Arguments
type(fv3jedi_geom)  , intent(in)  :: geom
real(kind=kind_real), intent(in)  :: day_of_year
real(kind=kind_real), intent(in)  :: field_slmsk          (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_sheleg         (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_tsea           (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_vtype          (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_stype          (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_vfrac          (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_stc            (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_smc            (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_u_srf          (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_v_srf          (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_f10m           (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_sss            (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: vegetation_type      (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: land_type_npoess     (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: land_type_igbp       (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: soil_type            (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: water_coverage       (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: land_coverage        (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: ice_coverage         (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: snow_coverage        (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: lai                  (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: water_temperature    (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: land_temperature     (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: ice_temperature      (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: snow_temperature     (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: soil_moisture_content(geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: vegetation_fraction  (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: soil_temperature     (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: snow_depth           (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: wind_speed           (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: wind_direction       (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: sea_surface_salinity (geom%isc:geom%iec,geom%jsc:geom%jec,1)

!Locals
real(kind=kind_real), parameter :: minsnow = 1.0_kind_real / 10.0_kind_real
real(kind=kind_real), parameter :: windlimit = 0.0001_kind_real
real(kind=kind_real), parameter :: quadcof  (4, 2  ) =      &
                                   reshape((/0.0_kind_real,  1.0_kind_real, 1.0_kind_real,  2.0_kind_real, &
                                             1.0_kind_real, -1.0_kind_real, 1.0_kind_real, -1.0_kind_real/), (/4, 2/))

integer              :: itype, istype
integer              :: istyp00
integer              :: lai_type, iquadrant
logical              :: lwind
real(kind=kind_real) :: sfcpct(0:3), ts(0:3), wgtavg(0:3), dtskin(0:3)
real(kind=kind_real) :: sno00
real(kind=kind_real) :: sst00
real(kind=kind_real) :: ss00
real(kind=kind_real) :: tsavg,ssavg
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
integer, parameter :: IGBP_N_TYPES = 20
 integer, parameter, dimension(1:IGBP_N_TYPES) :: igbp_to_npoess=(/PINE_FOREST, &
   BROADLEAF_FOREST, PINE_FOREST, BROADLEAF_FOREST, BROADLEAF_PINE_FOREST, &
   SCRUB, SCRUB_SOIL, BROADLEAF_BRUSH, BROADLEAF_BRUSH, SCRUB, BROADLEAF_BRUSH, &
   TILLED_SOIL, URBAN_CONCRETE, TILLED_SOIL, INVALID_LAND, COMPACTED_SOIL, &
   INVALID_LAND, TUNDRA, TUNDRA, TUNDRA/)
 integer, parameter, dimension(1:IGBP_N_TYPES) :: igbp_to_igbp=(/1, &
   2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, &
   20/)

 ! CRTM IR/vis uses 20 land surface types, but each type can be represented in any of three
 ! different classifications: NPOESS, IGBP, USGS. We use the GSI mappings from the GFS
 ! model type to IGBP and NPOESS. Currently the mapping to USGS is not implemented, but it could
 ! be added following exactly the same logic used in the other two cases.
 integer, parameter :: num_land_types = 20
 integer, parameter, dimension(1:num_land_types) :: map_model_sfc_to_crtm_land_npoess = &
   igbp_to_npoess
 integer, parameter, dimension(1:num_land_types) :: map_model_sfc_to_crtm_land_igbp = igbp_to_igbp

 ! CRTM microwave uses 13 vegetation types
 integer, parameter :: num_vege_types = 20
 integer, parameter, dimension(1:num_vege_types) :: map_model_sfc_to_crtm_mwave_vege=(/4, &
   1, 5, 2, 3, 8, 9, 6, 6, 7, 8, 12, 7, 12, 13, 11, 0, 10, 10, 11/)
 ! CRTM microwave uses 9 soil types
 integer, parameter :: num_soil_types = 16
 integer, parameter, dimension(1:num_soil_types) :: map_model_soil_to_crtm_mwave_soil=(/1, &
   1, 4, 2, 2, 8, 7, 2, 6, 5, 2, 3, 8, 1, 6, 9/)

integer :: ji, jj
integer              :: slmsk
integer              :: vtype
integer              :: stype
real(kind=kind_real) :: sheleg
real(kind=kind_real) :: tsea
real(kind=kind_real) :: vfrac
real(kind=kind_real) :: snwdph
real(kind=kind_real) :: stc
real(kind=kind_real) :: smc
real(kind=kind_real) :: u_srf
real(kind=kind_real) :: v_srf
real(kind=kind_real) :: f10m
real(kind=kind_real) :: sss

real(kind=kind_real) :: pi, rad2deg

! Constants
pi = constant('pi')
rad2deg = constant('rad2deg')

! Loop over all grid points
do jj = geom%jsc, geom%jec
  do ji = geom%isc, geom%iec

!   Why copy to scalars?
    slmsk  = nint(field_slmsk (ji,jj,1))
    vtype  = nint(field_vtype (ji,jj,1))
    stype  = nint(field_stype (ji,jj,1))
    sheleg = field_sheleg(ji,jj,1)
    tsea   = field_tsea  (ji,jj,1)
    vfrac  = field_vfrac (ji,jj,1)
    stc    = field_stc   (ji,jj,1)
    smc    = field_smc   (ji,jj,1)
    u_srf  = field_u_srf (ji,jj,1)
    v_srf  = field_v_srf (ji,jj,1)
    f10m   = field_f10m  (ji,jj,1)
    sss    = field_sss   (ji,jj,1)

    dtskin = 0.0_kind_real !TODO need real skin temperature increment?

    lwind = .true.

    ! Stage 1, like deter_sfc in GSI
    ! ------------------------------
    istyp00 = slmsk
    sno00 = sheleg !sno00 = snwdph
    sst00 = tsea
    tsavg = sst00
    ss00 = sss

    ssavg = ss00

    ! Potential for missing values in sno00 (if missing set to 0.0)
    if (abs(sno00) > 10.0e10_kind_real) sno00 = 0.0_kind_real

    if (istyp00 >=1 .and. sno00 > minsnow) istyp00 = 3

    sfcpct = 0.0_kind_real
    sfcpct(istyp00) = 1.0

    ts(0:3)=0.0_kind_real
    wgtavg(0:3)=0.0_kind_real
    vfr=0.0_kind_real
    stp=0.0_kind_real
    sty=0.0_kind_real
    vty=0.0_kind_real
    sm=0.0_kind_real
    sn=0.0_kind_real
    ss=0.0_kind_real

    if(istyp00 == 1)then
       vty  = vtype
       sty  = stype
       wgtavg(1) = 1.0
       ts(1)= sst00
       vfr  = vfrac
       stp  = stc
       sm   = smc
    else if(istyp00 == 2)then
       wgtavg(2) = 1.0
       ts(2)=sst00
    else if(istyp00 == 3)then
       wgtavg(3) = 1.0
       ts(3)=sst00
       sn = sno00
    else
       wgtavg(0) = 1.0
       ts(0)=sst00
       ss   =ss00
    end if

    if(wgtavg(0) > 0.0_kind_real)then
       ts(0) = ts(0)/wgtavg(0)
       ss    = ss/wgtavg(0)
    else
       ts(0) = tsavg
       ss    = ssavg
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

    f10 = f10m

    ! Stage 2 - like crtm_interface from GSI
    ! --------------------------------------
    ! If vty/sty will give maximum values, that will be out of range for CRTM, set to 1
    if (vty == 15) vty = 1
    if (sty == 16) sty = 1

    itype  = vty
    istype = sty

    itype  = min(max(1,itype),num_vege_types)
    istype = min(max(1,istype),num_soil_types)
    land_type_npoess(ji,jj,1) = real(max(1,map_model_sfc_to_crtm_land_npoess(itype)),kind_real)
    land_type_igbp(ji,jj,1) = real(max(1,map_model_sfc_to_crtm_land_igbp(itype)),kind_real)
    Vegetation_Type(ji,jj,1) = real(max(1,map_model_sfc_to_crtm_mwave_vege(itype)),kind_real)
    Soil_Type(ji,jj,1) = real(map_model_soil_to_crtm_mwave_soil(istype),kind_real)
    lai_type = real(map_model_sfc_to_crtm_mwave_vege(itype),kind_real)

    water_coverage(ji,jj,1) = min(max(0.0_kind_real,sfcpct(0)),1.0_kind_real)
    land_coverage(ji,jj,1)  = min(max(0.0_kind_real,sfcpct(1)),1.0_kind_real)
    ice_coverage(ji,jj,1)   = min(max(0.0_kind_real,sfcpct(2)),1.0_kind_real)
    snow_coverage(ji,jj,1)  = min(max(0.0_kind_real,sfcpct(3)),1.0_kind_real)

    lai(ji,jj,1) = 0.0_kind_real

    if (land_coverage(ji,jj,1) > 0.0_kind_real) then

       if (lai_type > 0) then
         call get_lai(lai_type, geom%grid_lat(ji,jj), day_of_year, lai(ji,jj,1))
       endif

       ! for Glacial land ice soil type and vegetation type
       if(Soil_Type(ji,jj,1) == 9 .OR. Vegetation_Type(ji,jj,1) == 13) then
          ice_coverage(ji,jj,1) = min(ice_coverage(ji,jj,1) + land_coverage(ji,jj,1), 1.0_kind_real)
          land_coverage(ji,jj,1) = 0.0_kind_real
       endif

    endif

    if (lwind) then

      !Interpolate lowest level winds to observation location
      uu5 = u_srf
      vv5 = v_srf
      windscale = 1.0_kind_real

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
      wind_speed(ji,jj,1)           = sfc_speed
      wind_direction(ji,jj,1)       = rad2deg*wind10_direction

    else

      wind_speed(ji,jj,1)           = 0.0_kind_real
      wind_direction(ji,jj,1)       = 0.0_kind_real

    endif

!   Why copy from scalars?
    water_temperature(ji,jj,1)     = max(ts(0) + dtskin(0), 270._kind_real)
    sea_surface_salinity(ji,jj,1)  = ss

    !TODO, is nst_gsi ever > 1?
    !if(nst_gsi > 1 .and. water_coverage(1) > 0.0_kind_real) then
       !water_temperature(ji,jj,1)  = max(data_s(itref)+data_s(idtw)-data_s(idtc) + dtskin(0), 271._kind_real)
    !endif

!   Why copy from scalars?
    land_temperature(ji,jj,1)      = ts(1) + dtskin(1)
    ice_temperature(ji,jj,1)       = min(ts(2) + dtskin(2), 280._kind_real)
    snow_temperature(ji,jj,1)      = min(ts(3) + dtskin(3), 280._kind_real)
    soil_moisture_content(ji,jj,1) = sm
    vegetation_fraction(ji,jj,1)   = vfr
    soil_temperature(ji,jj,1)      = stp
    snow_depth(ji,jj,1)            = sn

  enddo
enddo

end subroutine crtm_surface

!----------------------------------------------------------------------------

! Implement leaf-area index (LAI) from GSI's crtm_interface module
!
! This is a simple triangle-wave function reaching its min value in winter and its max value
! in summer (in each hemisphere). Note this means at the equator, where seasons are offset,
! the LAI is discontinuous, with min/winter values of LAI meeting max/summer values.
subroutine get_lai(lai_type, latitude, day_of_year_in, lai)

  implicit none

  integer, intent(in) :: lai_type
  real(kind=kind_real), intent(in) :: latitude
  real(kind=kind_real), intent(in) :: day_of_year_in
  real(kind=kind_real), intent(out) :: lai

  real(kind=kind_real), dimension(2) :: lai_season
  real(kind=kind_real), dimension(13), parameter :: lai_min = (/ &
          3.08_kind_real, 1.85_kind_real, 2.80_kind_real, 5.00_kind_real, 1.00_kind_real, &
          0.50_kind_real, 0.52_kind_real, 0.60_kind_real, 0.50_kind_real, 0.60_kind_real, &
          0.10_kind_real, 1.56_kind_real, 0.01_kind_real /)
  real(kind=kind_real), dimension(13), parameter :: lai_max = (/ &
          6.48_kind_real, 3.31_kind_real, 5.50_kind_real, 6.40_kind_real, 5.16_kind_real, &
          3.66_kind_real, 2.90_kind_real, 2.60_kind_real, 3.66_kind_real, 2.60_kind_real, &
          0.75_kind_real, 5.68_kind_real, 0.01_kind_real /)

  ! Days-of-year for mid-Jan & mid-Jul (leap years are ignored), when the simple LAI function
  ! reaches its extremal values
  real(kind=kind_real), dimension(3), parameter :: day_of_peak = &
          (/15.5_kind_real, 196.5_kind_real, 380.5_kind_real/)

  real(kind=kind_real) :: day_of_year
  integer :: ni, n1, n2
  real(kind=kind_real) :: w1, w2

  day_of_year = day_of_year_in
  if (day_of_year .lt. day_of_peak(1)) day_of_year = day_of_year + 365.0_kind_real

  do ni = 1,2
    if (day_of_year .ge. day_of_peak(ni) .and. day_of_year .lt. day_of_peak(ni + 1)) then
      n1 = ni
      n2 = ni + 1
      exit
    endif
    ! If got here, then day_of_year < day_of_peak(1) or day_of_year >= day_of_peak(3),
    ! so either day_of_year was very wrong on input or the rephasing logic above failed.
    if (ni == 2) then
      call abor1_ftn('fv3jedi.surface_vt_mod.get_lai received invalid day_of_year')
    endif
  enddo
  w1 = (day_of_peak(n2) - day_of_year) / (day_of_peak(n2) - day_of_peak(n1))
  w2 = (day_of_year - day_of_peak(n1)) / (day_of_peak(n2) - day_of_peak(n1))
  if (n2 .eq. 3) n2 = 1

  lai_season(1) = lai_min(lai_type)
  lai_season(2) = lai_max(lai_type)
  if (latitude < 0.0_kind_real) then
    lai = w1 * lai_season(n2) + w2 * lai_season(n1)
  else
    lai = w1 * lai_season(n1) + w2 * lai_season(n2)
  endif

end subroutine get_lai

!----------------------------------------------------------------------------

end module surface_vt_mod
