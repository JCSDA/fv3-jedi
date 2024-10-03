! (C) Copyright 2018-2019 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module surface_vt_mod

use fv3jedi_geom_mod, only: fv3jedi_geom
use fv3jedi_kinds_mod, only: kind_real
use fv3jedi_constants_mod, only: constant

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
real(kind=kind_real), parameter :: minswe = 1.0_kind_real / 10.0_kind_real
real(kind=kind_real), parameter :: windlimit = 0.0001_kind_real

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
integer, parameter :: num_vtypes = 20
integer, parameter, dimension(1:num_vtypes) :: map_model_sfc_to_crtm_land_npoess = &
  igbp_to_npoess
integer, parameter, dimension(1:num_vtypes) :: map_model_sfc_to_crtm_land_igbp = igbp_to_igbp

! CRTM microwave uses 13 vegetation types
integer, parameter, dimension(1:num_vtypes) :: map_model_sfc_to_crtm_mwave_vege=(/4, &
  1, 5, 2, 3, 8, 9, 6, 6, 7, 8, 12, 7, 12, 13, 11, 0, 10, 10, 11/)
! CRTM microwave uses 9 soil types
integer, parameter :: num_stypes = 16
integer, parameter, dimension(1:num_stypes) :: map_model_soil_to_crtm_mwave_soil=(/1, &
  1, 4, 2, 2, 8, 7, 2, 6, 5, 2, 3, 8, 1, 6, 9/)

real(kind=kind_real) :: local_swe(geom%isc:geom%iec,geom%jsc:geom%jec,1)
integer :: local_slmsk(geom%isc:geom%iec,geom%jsc:geom%jec,1)
integer :: ji, jj
integer :: vtype, stype, lai_veg_type
real(kind=kind_real) :: rad2deg

! Constants
rad2deg = constant('rad2deg')

! Potential for missing values in snow water equivalent (if missing set to 0.0)
local_swe = field_sheleg  ! SWE is named "sheleg" in backgrounds
where (abs(local_swe) > 10.0e10_kind_real) local_swe = 0.0_kind_real

! Redefine land/ice with snow => snow
! Note: The GFS slmsk has values {0,1,2} denoting {sea,land,ice}.
!       Locally within this function, we also use an additional value (3) to denote snow.
local_slmsk = nint(field_slmsk)
where (local_slmsk >= 1 .and. local_swe > minswe) local_slmsk = 3

! Defaults to override below
! Note that when requesting FOV averaging of surface fields from UFO CRTM, these coverage fractions
! will serve as the averaging masks for the different variables over each FOV. Care should be taken
! to ensure that the averaging masks passed to UFO CRTM are consistent with the interpolation masks
! passed to OOPS's GetValues for interpolation to obs locations. As a simple rule, the interpolation
! mask for a particular field should generally match the FOV averaging mask.
water_coverage = 0.0_kind_real
land_coverage = 0.0_kind_real
ice_coverage = 0.0_kind_real
snow_coverage = 0.0_kind_real

! When the interpolation to obs locations does NOT use surface-type masks, we have the unfortunate
! situation where surface-specific fields can be interpolated across surface types. This is a
! problem because there's no valid soil_type for ocean points (as an arbitrary example), so
! UNMASKED interpolation to an obs between land and ocean may produce an incorrect soil_type GeoVaL.
! The correct solution is to use surface-type interpolation masks to avoid this cross-contamination.
! Until masks are used by default in fv3-jedi experiments, we keep these GFS-based defaults for
! backward compatibility -- but not because they are somehow correct.
soil_temperature = 0.0_kind_real
soil_moisture_content = 1.0_kind_real
vegetation_fraction = 0.0_kind_real
land_type_npoess = 9.0_kind_real  ! pine forest
land_type_igbp = 1.0_kind_real  ! evergreen needleleaf forest
vegetation_type = 4.0_kind_real  ! evergreen needleleaf forest
soil_type = 1.0_kind_real  ! coarse loamy sand
lai = 0.0_kind_real
snow_depth = 0.0_kind_real

where (local_slmsk == 0)  ! water
  water_coverage = 1.0_kind_real
end where

where (local_slmsk == 1)  ! land
  land_coverage = 1.0_kind_real
  soil_temperature = field_stc
  soil_moisture_content = field_smc
  vegetation_fraction = field_vfrac
  ! veg type, land type, soil type, leaf area index are set in loops below
end where

where (local_slmsk == 2)  ! ice
  ice_coverage = 1.0_kind_real
end where

where (local_slmsk == 3)  ! local value indicates snow coverage (not in GFS backgrounds)
  snow_coverage = 1.0_kind_real
  snow_depth = local_swe  ! this assigns SWE to snow depth, probably an old GSI bug inherited here
end where

! These thresholds are (probably?) applied to the model fields to make them compatible with CRTM
! consistency checks on the temperatures for different surface types. A possible cleanup would be
! to move this check from fv3jedi to the UFO CRTM operator (before calling CRTM itself).
water_temperature = max(field_tsea, 270.0_kind_real)
land_temperature = field_tsea
ice_temperature = min(field_tsea, 280.0_kind_real)
snow_temperature = min(field_tsea, 280.0_kind_real)

wind_speed = field_f10m * sqrt(field_u_srf**2 + field_v_srf**2)
! atan2(y,x) gives rads north from east
! atan2(x,y) gives rads east from north, per CRTM definition
! convert to degrees and fix phasing to lie in [0,360]
wind_direction = rad2deg * atan2(field_u_srf, field_v_srf)
where (field_u_srf < 0.0_kind_real)
  wind_direction = wind_direction + 360.0_kind_real
end where

sea_surface_salinity = field_sss

! Loop over grid points
do jj = geom%jsc, geom%jec
  do ji = geom%isc, geom%iec
    if (local_slmsk(ji,jj,1) == 1) then

      vtype = nint(field_vtype(ji,jj,1))
      stype = nint(field_stype(ji,jj,1))

      ! For vtype / stype matching "glacial land ice" => reassign land to ice
      if (vtype == 15 .or. stype == 16) then
        ! set ice fields, reset land fields to default values
        ! (basically, set all these fields as if slmsk==2 from the beginning)
        ice_coverage(ji,jj,1) = 1.0_kind_real
        land_coverage(ji,jj,1) = 0.0_kind_real
        soil_temperature = 0.0_kind_real
        soil_moisture_content = 1.0_kind_real
        vegetation_fraction = 0.0_kind_real
        lai = 0.0_kind_real
        continue  ! work below assumes slmsk==1, so skip to next i,j
      end if

      ! This silently fixes out-of-range vtype/stype that would lead to indexing problems below.
      ! This means the code will always run, but if incorrect values are passed as vtype/stype, it
      ! will also be difficult to identify/debug the problem.
      vtype = min(max(1, vtype), num_vtypes)
      stype = min(max(1, stype), num_stypes)

      land_type_npoess(ji,jj,1) = real(max(1,map_model_sfc_to_crtm_land_npoess(vtype)), kind_real)
      land_type_igbp(ji,jj,1) = real(max(1,map_model_sfc_to_crtm_land_igbp(vtype)), kind_real)
      vegetation_type(ji,jj,1) = real(max(1,map_model_sfc_to_crtm_mwave_vege(vtype)), kind_real)
      soil_type(ji,jj,1) = real(map_model_soil_to_crtm_mwave_soil(stype), kind_real)

      lai_veg_type = nint(vegetation_type(ji,jj,1))
      if (lai_veg_type > 0) then
        call get_lai(lai_veg_type, geom%grid_lat(ji,jj), day_of_year, lai(ji,jj,1))
      endif
    endif  ! land vs sea/ice/snow
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
