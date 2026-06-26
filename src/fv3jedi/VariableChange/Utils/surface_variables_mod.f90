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
                         field_slmsk, field_sheleg, field_skin_temperature_at_surface, &
                         field_vtype, field_stype, field_vfrac, field_stc, field_smc, &
                         field_eastward_wind_at_surface, field_northward_wind_at_surface, &
                         field_wind_reduction_factor_at_10m, field_sea_surface_salinity, &
                         land_type_index_npoess, land_type_index_igbp, &
                         vegetation_type_index, soil_type, water_area_fraction, land_area_fraction, &
                         ice_area_fraction, surface_snow_area_fraction, leaf_area_index, &
                         skin_temperature_at_surface_where_sea, &
                         skin_temperature_at_surface_where_land, &
                         skin_temperature_at_surface_where_ice, &
                         skin_temperature_at_surface_where_snow, &
                         volume_fraction_of_condensed_water_in_soil, vegetation_area_fraction, &
                         soil_temperature, surface_snow_thickness, wind_speed_at_surface, &
                         wind_from_direction_at_surface, sea_surface_salinity)

implicit none

!Arguments
type(fv3jedi_geom)  , intent(in)  :: geom
real(kind=kind_real), intent(in)  :: day_of_year
real(kind=kind_real), intent(in)  :: field_slmsk                          (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_sheleg                         (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_skin_temperature_at_surface    (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_vtype                          (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_stype                          (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_vfrac                          (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_stc                            (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_smc                            (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_eastward_wind_at_surface       (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_northward_wind_at_surface      (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_wind_reduction_factor_at_10m   (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(in)  :: field_sea_surface_salinity           (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: vegetation_type_index                     (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: land_type_index_npoess                    (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: land_type_index_igbp                      (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: soil_type                                 (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: water_area_fraction                       (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: land_area_fraction                        (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: ice_area_fraction                         (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: surface_snow_area_fraction                (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: leaf_area_index                           (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: skin_temperature_at_surface_where_sea     (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: skin_temperature_at_surface_where_land    (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: skin_temperature_at_surface_where_ice     (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: skin_temperature_at_surface_where_snow    (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: volume_fraction_of_condensed_water_in_soil(geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: vegetation_area_fraction                  (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: soil_temperature                          (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: surface_snow_thickness                    (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: wind_speed_at_surface                     (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: wind_from_direction_at_surface            (geom%isc:geom%iec,geom%jsc:geom%jec,1)
real(kind=kind_real), intent(inout) :: sea_surface_salinity                      (geom%isc:geom%iec,geom%jsc:geom%jec,1)

! Parameters used by the crtm_surface_* routines below
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

!Locals
integer :: isc, iec, jsc, jec
real(kind=kind_real) :: local_swe  (geom%isc:geom%iec,geom%jsc:geom%jec,1)
integer              :: local_slmsk(geom%isc:geom%iec,geom%jsc:geom%jec,1)

isc = geom%isc
iec = geom%iec
jsc = geom%jsc
jec = geom%jec

! Shared intermediate quantities, each produced by its own subroutine.
call crtm_surface_local_swe(isc, iec, jsc, jec, field_sheleg, local_swe)
call crtm_surface_local_slmsk(isc, iec, jsc, jec, field_slmsk, local_swe, minswe, local_slmsk)

! Each call below produces exactly one output variable and never alters it again.
! Ordering constraint: vegetation_type_index must be produced before leaf_area_index, because
! leaf_area_index is derived from the final vegetation_type_index values.
call crtm_surface_water_coverage(isc, iec, jsc, jec, local_slmsk, water_area_fraction)
call crtm_surface_land_coverage(isc, iec, jsc, jec, local_slmsk, field_vtype, field_stype, &
                                land_area_fraction)
call crtm_surface_ice_coverage(isc, iec, jsc, jec, local_slmsk, field_vtype, field_stype, &
                               ice_area_fraction)
call crtm_surface_snow_coverage(isc, iec, jsc, jec, local_slmsk, surface_snow_area_fraction)
call crtm_surface_snow_depth(isc, iec, jsc, jec, local_slmsk, local_swe, surface_snow_thickness)
call crtm_surface_soil_temperature(isc, iec, jsc, jec, local_slmsk, field_stc, field_vtype, &
                                   field_stype, soil_temperature)
call crtm_surface_soil_moisture_content(isc, iec, jsc, jec, local_slmsk, field_smc, field_vtype, &
                                        field_stype, volume_fraction_of_condensed_water_in_soil)
call crtm_surface_vegetation_fraction(isc, iec, jsc, jec, local_slmsk, field_vfrac, field_vtype, &
                                      field_stype, vegetation_area_fraction)
call crtm_surface_water_temperature(isc, iec, jsc, jec, field_skin_temperature_at_surface, &
                                    skin_temperature_at_surface_where_sea)
call crtm_surface_land_temperature(isc, iec, jsc, jec, field_skin_temperature_at_surface, &
                                   skin_temperature_at_surface_where_land)
call crtm_surface_ice_temperature(isc, iec, jsc, jec, field_skin_temperature_at_surface, &
                                  skin_temperature_at_surface_where_ice)
call crtm_surface_snow_temperature(isc, iec, jsc, jec, field_skin_temperature_at_surface, &
                                   skin_temperature_at_surface_where_snow)
call crtm_surface_wind_speed(isc, iec, jsc, jec, field_wind_reduction_factor_at_10m, &
                             field_eastward_wind_at_surface, field_northward_wind_at_surface, &
                             wind_speed_at_surface)
call crtm_surface_wind_direction(isc, iec, jsc, jec, field_eastward_wind_at_surface, &
                                 field_northward_wind_at_surface, wind_from_direction_at_surface)
call crtm_surface_sea_surface_salinity(isc, iec, jsc, jec, field_sea_surface_salinity, &
                                       sea_surface_salinity)
call crtm_surface_land_type_npoess(isc, iec, jsc, jec, local_slmsk, field_vtype, &
                                   map_model_sfc_to_crtm_land_npoess, land_type_index_npoess)
call crtm_surface_land_type_igbp(isc, iec, jsc, jec, local_slmsk, field_vtype, &
                                 map_model_sfc_to_crtm_land_igbp, land_type_index_igbp)
call crtm_surface_vegetation_type(isc, iec, jsc, jec, local_slmsk, field_vtype, &
                                  map_model_sfc_to_crtm_mwave_vege, vegetation_type_index)
call crtm_surface_soil_type(isc, iec, jsc, jec, local_slmsk, field_stype, &
                            map_model_soil_to_crtm_mwave_soil, soil_type)
call crtm_surface_lai(isc, iec, jsc, jec, local_slmsk, field_vtype, field_stype, &
                      vegetation_type_index, geom%grid_lat(isc:iec,jsc:jec), day_of_year, &
                      leaf_area_index)

end subroutine crtm_surface

!----------------------------------------------------------------------------
! Shared intermediates -------------------------------------------------------
!----------------------------------------------------------------------------

! Snow water equivalent, with missing values (if any) reset to 0.0
subroutine crtm_surface_local_swe(isc, iec, jsc, jec, field_sheleg, local_swe)
integer             , intent(in)  :: isc, iec, jsc, jec
real(kind=kind_real), intent(in)  :: field_sheleg(isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: local_swe   (isc:iec,jsc:jec,1)

local_swe = field_sheleg  ! SWE is named "sheleg" in backgrounds
where (abs(local_swe) > 10.0e10_kind_real) local_swe = 0.0_kind_real

end subroutine crtm_surface_local_swe

!----------------------------------------------------------------------------

! Redefine land/ice with snow => snow
! Note: The GFS slmsk has values {0,1,2} denoting {sea,land,ice}.
!       Locally within these routines, we also use an additional value (3) to denote snow.
subroutine crtm_surface_local_slmsk(isc, iec, jsc, jec, field_slmsk, local_swe, minswe, &
                                    local_slmsk)
integer             , intent(in)  :: isc, iec, jsc, jec
real(kind=kind_real), intent(in)  :: field_slmsk(isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: local_swe  (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: minswe
integer             , intent(out) :: local_slmsk(isc:iec,jsc:jec,1)

local_slmsk = nint(field_slmsk)
where (local_slmsk >= 1 .and. local_swe > minswe) local_slmsk = 3

end subroutine crtm_surface_local_slmsk

!----------------------------------------------------------------------------
! Surface coverage fractions -------------------------------------------------
!----------------------------------------------------------------------------

subroutine crtm_surface_water_coverage(isc, iec, jsc, jec, local_slmsk, water_area_fraction)
integer             , intent(in)  :: isc, iec, jsc, jec
integer             , intent(in)  :: local_slmsk        (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: water_area_fraction(isc:iec,jsc:jec,1)

water_area_fraction = 0.0_kind_real
where (local_slmsk == 0) water_area_fraction = 1.0_kind_real

end subroutine crtm_surface_water_coverage

!----------------------------------------------------------------------------

subroutine crtm_surface_land_coverage(isc, iec, jsc, jec, local_slmsk, field_vtype, field_stype, &
                                      land_area_fraction)
integer             , intent(in)  :: isc, iec, jsc, jec
integer             , intent(in)  :: local_slmsk       (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_vtype       (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_stype       (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: land_area_fraction(isc:iec,jsc:jec,1)
integer :: ji, jj, vtype, stype

land_area_fraction = 0.0_kind_real
where (local_slmsk == 1) land_area_fraction = 1.0_kind_real

! For vtype / stype matching "glacial land ice" => reassign land to ice
do jj = jsc, jec
  do ji = isc, iec
    if (local_slmsk(ji,jj,1) == 1) then
      vtype = nint(field_vtype(ji,jj,1))
      stype = nint(field_stype(ji,jj,1))
      if (vtype == 15 .or. stype == 16) land_area_fraction(ji,jj,1) = 0.0_kind_real
    endif
  enddo
enddo

end subroutine crtm_surface_land_coverage

!----------------------------------------------------------------------------

subroutine crtm_surface_ice_coverage(isc, iec, jsc, jec, local_slmsk, field_vtype, field_stype, &
                                     ice_area_fraction)
integer             , intent(in)  :: isc, iec, jsc, jec
integer             , intent(in)  :: local_slmsk      (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_vtype      (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_stype      (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: ice_area_fraction(isc:iec,jsc:jec,1)
integer :: ji, jj, vtype, stype

ice_area_fraction = 0.0_kind_real
where (local_slmsk == 2) ice_area_fraction = 1.0_kind_real

! For vtype / stype matching "glacial land ice" => reassign land to ice
do jj = jsc, jec
  do ji = isc, iec
    if (local_slmsk(ji,jj,1) == 1) then
      vtype = nint(field_vtype(ji,jj,1))
      stype = nint(field_stype(ji,jj,1))
      if (vtype == 15 .or. stype == 16) ice_area_fraction(ji,jj,1) = 1.0_kind_real
    endif
  enddo
enddo

end subroutine crtm_surface_ice_coverage

!----------------------------------------------------------------------------

subroutine crtm_surface_snow_coverage(isc, iec, jsc, jec, local_slmsk, surface_snow_area_fraction)
integer             , intent(in)  :: isc, iec, jsc, jec
integer             , intent(in)  :: local_slmsk               (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: surface_snow_area_fraction(isc:iec,jsc:jec,1)

surface_snow_area_fraction = 0.0_kind_real
! local slmsk value of 3 indicates snow coverage
where (local_slmsk == 3) surface_snow_area_fraction = 1.0_kind_real

end subroutine crtm_surface_snow_coverage

!----------------------------------------------------------------------------
! Surface fields -------------------------------------------------------------
!----------------------------------------------------------------------------

subroutine crtm_surface_snow_depth(isc, iec, jsc, jec, local_slmsk, local_swe, &
                                   surface_snow_thickness)
integer             , intent(in)  :: isc, iec, jsc, jec
integer             , intent(in)  :: local_slmsk           (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: local_swe             (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: surface_snow_thickness(isc:iec,jsc:jec,1)

surface_snow_thickness = 0.0_kind_real
! this assigns SWE to snow depth, probably an old GSI bug inherited here
where (local_slmsk == 3) surface_snow_thickness = local_swe

end subroutine crtm_surface_snow_depth

!----------------------------------------------------------------------------

subroutine crtm_surface_soil_temperature(isc, iec, jsc, jec, local_slmsk, field_stc, &
                                         field_vtype, field_stype, soil_temperature)
integer             , intent(in)  :: isc, iec, jsc, jec
integer             , intent(in)  :: local_slmsk     (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_stc       (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_vtype     (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_stype     (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: soil_temperature(isc:iec,jsc:jec,1)
integer :: ji, jj, vtype, stype

soil_temperature = 0.0_kind_real
where (local_slmsk == 1) soil_temperature = field_stc

! For vtype / stype matching "glacial land ice" => reset to default value
do jj = jsc, jec
  do ji = isc, iec
    if (local_slmsk(ji,jj,1) == 1) then
      vtype = nint(field_vtype(ji,jj,1))
      stype = nint(field_stype(ji,jj,1))
      if (vtype == 15 .or. stype == 16) soil_temperature = 0.0_kind_real
    endif
  enddo
enddo

end subroutine crtm_surface_soil_temperature

!----------------------------------------------------------------------------

subroutine crtm_surface_soil_moisture_content(isc, iec, jsc, jec, local_slmsk, field_smc, &
                                              field_vtype, field_stype, &
                                              volume_fraction_of_condensed_water_in_soil)
integer             , intent(in)  :: isc, iec, jsc, jec
integer             , intent(in)  :: local_slmsk                               (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_smc                                 (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_vtype                               (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_stype                               (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: volume_fraction_of_condensed_water_in_soil(isc:iec,jsc:jec,1)
integer :: ji, jj, vtype, stype

volume_fraction_of_condensed_water_in_soil = 1.0_kind_real
where (local_slmsk == 1) volume_fraction_of_condensed_water_in_soil = field_smc

! For vtype / stype matching "glacial land ice" => reset to default value
do jj = jsc, jec
  do ji = isc, iec
    if (local_slmsk(ji,jj,1) == 1) then
      vtype = nint(field_vtype(ji,jj,1))
      stype = nint(field_stype(ji,jj,1))
      if (vtype == 15 .or. stype == 16) volume_fraction_of_condensed_water_in_soil = 1.0_kind_real
    endif
  enddo
enddo

end subroutine crtm_surface_soil_moisture_content

!----------------------------------------------------------------------------

subroutine crtm_surface_vegetation_fraction(isc, iec, jsc, jec, local_slmsk, field_vfrac, &
                                            field_vtype, field_stype, vegetation_area_fraction)
integer             , intent(in)  :: isc, iec, jsc, jec
integer             , intent(in)  :: local_slmsk             (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_vfrac             (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_vtype             (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_stype             (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: vegetation_area_fraction(isc:iec,jsc:jec,1)
integer :: ji, jj, vtype, stype

vegetation_area_fraction = 0.0_kind_real
where (local_slmsk == 1) vegetation_area_fraction = field_vfrac

! For vtype / stype matching "glacial land ice" => reset to default value
do jj = jsc, jec
  do ji = isc, iec
    if (local_slmsk(ji,jj,1) == 1) then
      vtype = nint(field_vtype(ji,jj,1))
      stype = nint(field_stype(ji,jj,1))
      if (vtype == 15 .or. stype == 16) vegetation_area_fraction = 0.0_kind_real
    endif
  enddo
enddo

end subroutine crtm_surface_vegetation_fraction

!----------------------------------------------------------------------------
! Surface temperatures -------------------------------------------------------
! These thresholds are (probably?) applied to the model fields to make them compatible with CRTM
! consistency checks on the temperatures for different surface types.
!----------------------------------------------------------------------------

subroutine crtm_surface_water_temperature(isc, iec, jsc, jec, field_skin_temperature_at_surface, &
                                          skin_temperature_at_surface_where_sea)
integer             , intent(in)  :: isc, iec, jsc, jec
real(kind=kind_real), intent(in)  :: field_skin_temperature_at_surface    (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: skin_temperature_at_surface_where_sea(isc:iec,jsc:jec,1)

skin_temperature_at_surface_where_sea = max(field_skin_temperature_at_surface, 270.0_kind_real)

end subroutine crtm_surface_water_temperature

!----------------------------------------------------------------------------

subroutine crtm_surface_land_temperature(isc, iec, jsc, jec, field_skin_temperature_at_surface, &
                                         skin_temperature_at_surface_where_land)
integer             , intent(in)  :: isc, iec, jsc, jec
real(kind=kind_real), intent(in)  :: field_skin_temperature_at_surface     (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: skin_temperature_at_surface_where_land(isc:iec,jsc:jec,1)

skin_temperature_at_surface_where_land = field_skin_temperature_at_surface

end subroutine crtm_surface_land_temperature

!----------------------------------------------------------------------------

subroutine crtm_surface_ice_temperature(isc, iec, jsc, jec, field_skin_temperature_at_surface, &
                                        skin_temperature_at_surface_where_ice)
integer             , intent(in)  :: isc, iec, jsc, jec
real(kind=kind_real), intent(in)  :: field_skin_temperature_at_surface    (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: skin_temperature_at_surface_where_ice(isc:iec,jsc:jec,1)

skin_temperature_at_surface_where_ice = min(field_skin_temperature_at_surface, 280.0_kind_real)

end subroutine crtm_surface_ice_temperature

!----------------------------------------------------------------------------

subroutine crtm_surface_snow_temperature(isc, iec, jsc, jec, field_skin_temperature_at_surface, &
                                         skin_temperature_at_surface_where_snow)
integer             , intent(in)  :: isc, iec, jsc, jec
real(kind=kind_real), intent(in)  :: field_skin_temperature_at_surface     (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: skin_temperature_at_surface_where_snow(isc:iec,jsc:jec,1)

skin_temperature_at_surface_where_snow = min(field_skin_temperature_at_surface, 280.0_kind_real)

end subroutine crtm_surface_snow_temperature

!----------------------------------------------------------------------------
! Surface winds and salinity -------------------------------------------------
!----------------------------------------------------------------------------

subroutine crtm_surface_wind_speed(isc, iec, jsc, jec, field_wind_reduction_factor_at_10m, &
                                   field_eastward_wind_at_surface, &
                                   field_northward_wind_at_surface, wind_speed_at_surface)
integer             , intent(in)  :: isc, iec, jsc, jec
real(kind=kind_real), intent(in)  :: field_wind_reduction_factor_at_10m(isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_eastward_wind_at_surface    (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_northward_wind_at_surface   (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: wind_speed_at_surface             (isc:iec,jsc:jec,1)

wind_speed_at_surface = field_wind_reduction_factor_at_10m * &
     sqrt(field_eastward_wind_at_surface**2 + field_northward_wind_at_surface**2)

end subroutine crtm_surface_wind_speed

!----------------------------------------------------------------------------

subroutine crtm_surface_wind_direction(isc, iec, jsc, jec, field_eastward_wind_at_surface, &
                                       field_northward_wind_at_surface, &
                                       wind_from_direction_at_surface)
integer             , intent(in)  :: isc, iec, jsc, jec
real(kind=kind_real), intent(in)  :: field_eastward_wind_at_surface (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_northward_wind_at_surface(isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: wind_from_direction_at_surface (isc:iec,jsc:jec,1)
real(kind=kind_real) :: rad2deg

rad2deg = constant('rad2deg')

! atan2(y,x) gives rads north from east
! atan2(x,y) gives rads east from north, per CRTM definition
! convert to degrees and fix phasing to lie in [0,360]
wind_from_direction_at_surface = rad2deg * atan2(field_eastward_wind_at_surface, &
                                                 field_northward_wind_at_surface)
where (field_eastward_wind_at_surface < 0.0_kind_real)
  wind_from_direction_at_surface = wind_from_direction_at_surface + 360.0_kind_real
end where

end subroutine crtm_surface_wind_direction

!----------------------------------------------------------------------------

subroutine crtm_surface_sea_surface_salinity(isc, iec, jsc, jec, field_sea_surface_salinity, &
                                             sea_surface_salinity)
integer             , intent(in)  :: isc, iec, jsc, jec
real(kind=kind_real), intent(in)  :: field_sea_surface_salinity(isc:iec,jsc:jec,1)
real(kind=kind_real), intent(out) :: sea_surface_salinity      (isc:iec,jsc:jec,1)

sea_surface_salinity = field_sea_surface_salinity

end subroutine crtm_surface_sea_surface_salinity

!----------------------------------------------------------------------------
! Surface / vegetation / soil classification ---------------------------------
!----------------------------------------------------------------------------

subroutine crtm_surface_land_type_npoess(isc, iec, jsc, jec, local_slmsk, field_vtype, &
                                         map_model_sfc_to_crtm_land_npoess, land_type_index_npoess)
integer             , intent(in)  :: isc, iec, jsc, jec
integer             , intent(in)  :: local_slmsk           (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_vtype           (isc:iec,jsc:jec,1)
integer             , intent(in)  :: map_model_sfc_to_crtm_land_npoess(:)
real(kind=kind_real), intent(out) :: land_type_index_npoess(isc:iec,jsc:jec,1)
integer :: ji, jj, vtype

land_type_index_npoess = 9.0_kind_real  ! pine forest

do jj = jsc, jec
  do ji = isc, iec
    if (local_slmsk(ji,jj,1) == 1) then
      ! Silently fix out-of-range vtype that would lead to indexing problems below.
      vtype = min(max(1, nint(field_vtype(ji,jj,1))), size(map_model_sfc_to_crtm_land_npoess))
      land_type_index_npoess(ji,jj,1) = &
           real(max(1,map_model_sfc_to_crtm_land_npoess(vtype)), kind_real)
    endif
  enddo
enddo

end subroutine crtm_surface_land_type_npoess

!----------------------------------------------------------------------------

subroutine crtm_surface_land_type_igbp(isc, iec, jsc, jec, local_slmsk, field_vtype, &
                                       map_model_sfc_to_crtm_land_igbp, land_type_index_igbp)
integer             , intent(in)  :: isc, iec, jsc, jec
integer             , intent(in)  :: local_slmsk         (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_vtype         (isc:iec,jsc:jec,1)
integer             , intent(in)  :: map_model_sfc_to_crtm_land_igbp(:)
real(kind=kind_real), intent(out) :: land_type_index_igbp(isc:iec,jsc:jec,1)
integer :: ji, jj, vtype

land_type_index_igbp = 1.0_kind_real  ! evergreen needleleaf forest

do jj = jsc, jec
  do ji = isc, iec
    if (local_slmsk(ji,jj,1) == 1) then
      ! Silently fix out-of-range vtype that would lead to indexing problems below.
      vtype = min(max(1, nint(field_vtype(ji,jj,1))), size(map_model_sfc_to_crtm_land_igbp))
      land_type_index_igbp(ji,jj,1) = &
           real(max(1,map_model_sfc_to_crtm_land_igbp(vtype)), kind_real)
    endif
  enddo
enddo

end subroutine crtm_surface_land_type_igbp

!----------------------------------------------------------------------------

subroutine crtm_surface_vegetation_type(isc, iec, jsc, jec, local_slmsk, field_vtype, &
                                        map_model_sfc_to_crtm_mwave_vege, vegetation_type_index)
integer             , intent(in)  :: isc, iec, jsc, jec
integer             , intent(in)  :: local_slmsk          (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_vtype          (isc:iec,jsc:jec,1)
integer             , intent(in)  :: map_model_sfc_to_crtm_mwave_vege(:)
real(kind=kind_real), intent(out) :: vegetation_type_index(isc:iec,jsc:jec,1)
integer :: ji, jj, vtype

vegetation_type_index = 4.0_kind_real  ! evergreen needleleaf forest

do jj = jsc, jec
  do ji = isc, iec
    if (local_slmsk(ji,jj,1) == 1) then
      ! Silently fix out-of-range vtype that would lead to indexing problems below.
      vtype = min(max(1, nint(field_vtype(ji,jj,1))), size(map_model_sfc_to_crtm_mwave_vege))
      vegetation_type_index(ji,jj,1) = &
           real(max(1,map_model_sfc_to_crtm_mwave_vege(vtype)), kind_real)
    endif
  enddo
enddo

end subroutine crtm_surface_vegetation_type

!----------------------------------------------------------------------------

subroutine crtm_surface_soil_type(isc, iec, jsc, jec, local_slmsk, field_stype, &
                                  map_model_soil_to_crtm_mwave_soil, soil_type)
integer             , intent(in)  :: isc, iec, jsc, jec
integer             , intent(in)  :: local_slmsk(isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_stype(isc:iec,jsc:jec,1)
integer             , intent(in)  :: map_model_soil_to_crtm_mwave_soil(:)
real(kind=kind_real), intent(out) :: soil_type  (isc:iec,jsc:jec,1)
integer :: ji, jj, stype

soil_type = 1.0_kind_real  ! coarse loamy sand

do jj = jsc, jec
  do ji = isc, iec
    if (local_slmsk(ji,jj,1) == 1) then
      ! Silently fix out-of-range stype that would lead to indexing problems below.
      stype = min(max(1, nint(field_stype(ji,jj,1))), size(map_model_soil_to_crtm_mwave_soil))
      soil_type(ji,jj,1) = real(map_model_soil_to_crtm_mwave_soil(stype), kind_real)
    endif
  enddo
enddo

end subroutine crtm_surface_soil_type

!----------------------------------------------------------------------------

! Leaf area index. Derived from the final vegetation_type_index, so vegetation_type_index must
! be produced by crtm_surface_vegetation_type before this routine is called.
!
! Implements leaf-area index (LAI) from GSI's crtm_interface module: a simple triangle-wave
! function reaching its min value in winter and its max value in summer (in each hemisphere).
! Note this means at the equator, where seasons are offset, the LAI is discontinuous, with
! min/winter values of LAI meeting max/summer values.
subroutine crtm_surface_lai(isc, iec, jsc, jec, local_slmsk, field_vtype, field_stype, &
                            vegetation_type_index, grid_lat, day_of_year, leaf_area_index)
integer             , intent(in)  :: isc, iec, jsc, jec
integer             , intent(in)  :: local_slmsk          (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_vtype          (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: field_stype          (isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: vegetation_type_index(isc:iec,jsc:jec,1)
real(kind=kind_real), intent(in)  :: grid_lat             (isc:iec,jsc:jec)
real(kind=kind_real), intent(in)  :: day_of_year
real(kind=kind_real), intent(out) :: leaf_area_index      (isc:iec,jsc:jec,1)

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

integer              :: ji, jj, vtype, stype, lai_veg_type
integer              :: ni, n1, n2
real(kind=kind_real) :: doy, w1, w2
real(kind=kind_real) :: lai_season(2)

leaf_area_index = 0.0_kind_real

! Compute the season interpolation weights once (they do not depend on the cell).
doy = day_of_year
if (doy .lt. day_of_peak(1)) doy = doy + 365.0_kind_real

n1 = 0
do ni = 1,2
  if (doy .ge. day_of_peak(ni) .and. doy .lt. day_of_peak(ni + 1)) then
    n1 = ni
    n2 = ni + 1
    exit
  endif
enddo
! If n1 is still 0, then doy < day_of_peak(1) or doy >= day_of_peak(3), so either day_of_year
! was very wrong on input or the rephasing logic above failed.
if (n1 == 0) then
  call abor1_ftn('fv3jedi.surface_vt_mod.crtm_surface_lai received invalid day_of_year')
endif
w1 = (day_of_peak(n2) - doy) / (day_of_peak(n2) - day_of_peak(n1))
w2 = (doy - day_of_peak(n1)) / (day_of_peak(n2) - day_of_peak(n1))
if (n2 .eq. 3) n2 = 1

do jj = jsc, jec
  do ji = isc, iec
    if (local_slmsk(ji,jj,1) == 1) then
      vtype = nint(field_vtype(ji,jj,1))
      stype = nint(field_stype(ji,jj,1))
      ! For vtype / stype matching "glacial land ice" => reset to default value
      if (vtype == 15 .or. stype == 16) leaf_area_index = 0.0_kind_real
      lai_veg_type = nint(vegetation_type_index(ji,jj,1))
      if (lai_veg_type > 0) then
        lai_season(1) = lai_min(lai_veg_type)
        lai_season(2) = lai_max(lai_veg_type)
        if (grid_lat(ji,jj) < 0.0_kind_real) then
          leaf_area_index(ji,jj,1) = w1 * lai_season(n2) + w2 * lai_season(n1)
        else
          leaf_area_index(ji,jj,1) = w1 * lai_season(n1) + w2 * lai_season(n2)
        endif
      endif
    endif
  enddo
enddo

end subroutine crtm_surface_lai

!----------------------------------------------------------------------------

end module surface_vt_mod
