/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#pragma once

#include <map>
#include <string>
#include <utility>

#include "fv3jedi/FieldMetadata/FieldsMetadata.h"

namespace fv3jedi {

  // -----------------------------------------------------------------------------------------------

  // Elements to be populated for each field
  struct metadataStruct {
    std::string longName;
    std::string shortName;
    std::string units;
    std::string kind;
    std::string tracer;  // Turned into bool but input as string to allow default checking
    std::string horizontalStaggerLocation;
    std::string levels;
    std::string space;
  };

  // -----------------------------------------------------------------------------------------------

  void setStructDefaults(struct metadataStruct & md) {
    md.longName = "long name";
    md.shortName = "short name";
    md.units = "units";
    md.kind = "kind";
    md.tracer = "tracer";
    md.horizontalStaggerLocation = "horizontal stagger location";
    md.levels = "levels";
    md.space = "space";
  }

  // -----------------------------------------------------------------------------------------------

  void assertStructNotDefaults(struct metadataStruct & md) {
    // Check that structure does not contain defaults
    // ----------------------------------------------
    ASSERT_MSG(md.longName != "long name", "long name was not set");
    ASSERT_MSG(md.shortName != "short name", "short name was not set");
    ASSERT_MSG(md.units != "units", "units was not set");
    ASSERT_MSG(md.kind != "kind", "kind was not set");
    ASSERT_MSG(md.tracer != "tracer", "tracer was not set");
    ASSERT_MSG(md.horizontalStaggerLocation != "horizontal stagger location", "stagger not set");
    ASSERT_MSG(md.levels != "levels", "levels was not set");
    ASSERT_MSG(md.space != "space", "space was not set");
  }

  // -----------------------------------------------------------------------------------------------

  void addFieldMetadata(std::map<std::string, FieldMetadata> & fieldsmetadata, const int & nlev,
                        struct metadataStruct & md) {
    // Check that structure does not contain defaults
    assertStructNotDefaults(md);

    // Create object to hold the metadata for this field
    FieldMetadata fieldmetadata(md.longName, nlev);

    // Populate the object
    fieldmetadata.setShrtName(md.shortName);
    fieldmetadata.setVarUnits(md.units);
    fieldmetadata.setDataKind(md.kind);
    fieldmetadata.setStagrLoc(md.horizontalStaggerLocation);
    fieldmetadata.setNumLevls(md.levels);
    fieldmetadata.setMathSpac(md.space);
    fieldmetadata.setIsTracer(md.tracer);
    fieldmetadata.setInOuName(md.shortName);  // Default to short, can only be set using override
    fieldmetadata.setInOuFile("default");     // Can only be set using override
    fieldmetadata.setIntrpTyp("default");     // Can only be set using override
    fieldmetadata.setIntrpMsk("default");     // Can only be set using override

    // Validate the choices
    fieldmetadata.validate();

    // Check key not already in the map
    ASSERT_MSG(fieldsmetadata.find(md.longName) == fieldsmetadata.end(),
               "FieldMetadataDefault::addFieldMetadata: Long name "+md.longName+" already used.");

    // Insert the object into the map
    fieldsmetadata.insert(std::pair<std::string, FieldMetadata>(md.longName, fieldmetadata));

    // Set back to defaults
    setStructDefaults(md);
  }

  // -----------------------------------------------------------------------------------------------

  void setDefaults(std::map<std::string, FieldMetadata> & fieldsmetadata, const int & nlev) {
    // Create structure and set to default
    struct metadataStruct md;
    setStructDefaults(md);

    // Default field metadata
    // ----------------------

    md.longName = "u_component_of_native_D_grid_wind";
    md.shortName = "ud";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "northsouth";
    md.levels = "full";
    md.space = "vector";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "v_component_of_native_D_grid_wind";
    md.shortName = "vd";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "eastwest";
    md.levels = "full";
    md.space = "vector";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "u_component_of_native_C_grid_wind";
    md.shortName = "uc";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "eastwest";
    md.levels = "full";
    md.space = "vector";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "v_component_of_native_C_grid_wind";
    md.shortName = "vc";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "northsouth";
    md.levels = "full";
    md.space = "vector";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "eastward_wind";
    md.shortName = "ua";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "northward_wind";
    md.shortName = "va";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "air_temperature";
    md.shortName = "t";
    md.units = "K";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "virtual_temperature";
    md.shortName = "tv";
    md.units = "K";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "air_potential_temperature";
    md.shortName = "pt";
    md.units = "K";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "air_pressure_thickness";
    md.shortName = "delp";
    md.units = "pa";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "air_pressure_to_kappa";
    md.shortName = "pkz";
    md.units = "Pa";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "air_pressure_levels";
    md.shortName = "pe";
    md.units = "Pa";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "air_pressure";
    md.shortName = "p";
    md.units = "Pa";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "air_pressure_at_surface";
    md.shortName = "ps";
    md.units = "Pa";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "ln_air_pressure_at_interface";
    md.shortName = "lnpe";
    md.units = "Pa";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "water_vapor_mixing_ratio_wrt_moist_air";
    md.shortName = "sphum";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "saturation_specific_humidity";
    md.shortName = "qsat";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "relative_humidity";
    md.shortName = "rh";
    md.units = "1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "cloud_liquid_ice";
    md.shortName = "ice_wat";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "cloud_liquid_water";
    md.shortName = "liq_wat";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_large_scale_cloud_ice_water";
    md.shortName = "qils";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_convective_cloud_ice_water";
    md.shortName = "qicn";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_large_scale_cloud_liquid_water";
    md.shortName = "qlls";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_convective_cloud_liquid_water";
    md.shortName = "qlcn";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "fraction_of_large_scale_cloud_that_is_ice";
    md.shortName = "qilsf";
    md.units = "1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "fraction_of_convective_cloud_that_is_ice";
    md.shortName = "qicnf";
    md.units = "1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "snow_water";
    md.shortName = "snowwat";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "rain_water";
    md.shortName = "rainwat";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "graupel";
    md.shortName = "graupel";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "cloud_droplet_number_concentration";
    md.shortName = "water_nc";
    md.units = "kg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "cloud_ice_number_concentration";
    md.shortName = "ice_nc";
    md.units = "kg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "rain_number_concentration";
    md.shortName = "rain_nc";
    md.units = "kg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "aerosol_water_number_concentration";
    md.shortName = "liq_aero";
    md.units = "kg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "aerosol_ice_number_concentration";
    md.shortName = "ice_aero";
    md.units = "kg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "cloud_area_fraction_in_atmosphere_layer";
    md.shortName = "cld_amt";
    md.units = "1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "sgs_tke";
    md.shortName = "sgs_tke";
    md.units = "m2/s2";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "upward_air_velocity";
    md.shortName = "w";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "layer_thickness";
    md.shortName = "delz";
    md.units = "m";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "geopotential_height_times_gravity_at_surface";
    md.shortName = "phis";
    md.units = "m";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "stream_function";
    md.shortName = "psi";
    md.units = "m+2s";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "velocity_potential";
    md.shortName = "chi";
    md.units = "m+2s";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "vorticity";
    md.shortName = "vort";
    md.units = "m+2s";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "divergence";
    md.shortName = "divg";
    md.units = "m+2s";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "slmsk";
    md.shortName = "slmsk";
    md.units = "none";
    md.kind = "integer";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "sheleg";
    md.shortName = "sheleg";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "skin_temperature_at_surface";
    md.shortName = "ts";
    md.units = "K";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "sea_surface_temperature";
    md.shortName = "sst";
    md.units = "K";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "vtype";
    md.shortName = "vtype";
    md.units = "none";
    md.kind = "integer";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "stype";
    md.shortName = "stype";
    md.units = "none";
    md.kind = "integer";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "vfrac";
    md.shortName = "vfrac";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "stc";
    md.shortName = "stc";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "4";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "tslb";
    md.shortName = "tslb";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "9";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "soilt";
    md.shortName = "soilt";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "soilMoistureVolumetric";
    md.shortName = "smc";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "4";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "smois";
    md.shortName = "smois";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "9";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "soilm";
    md.shortName = "soilm";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "totalSnowDepth";
    md.shortName = "snwdph";
    md.units = "mm";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "totalSnowDepthMeters";
    md.shortName = "snwdphMeters";
    md.units = "m";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "eastward_wind_at_surface";
    md.shortName = "u_srf";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "northward_wind_at_surface";
    md.shortName = "v_srf";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "f10m";
    md.shortName = "f10m";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "sea_surface_salinity";
    md.shortName = "sss";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "initial_mass_fraction_of_large_scale_cloud_condensate";
    md.shortName = "qls";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "initial_mass_fraction_of_convective_cloud_condensate";
    md.shortName = "qcn";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "convective_cloud_area_fraction";
    md.shortName = "cfcn";
    md.units = "1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "fraction_of_ocean";
    md.shortName = "frocean";
    md.units = "1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "fraction_of_land";
    md.shortName = "frland";
    md.units = "1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "fraction_of_landice";
    md.shortName = "frlandice";
    md.units = "1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "fraction_of_lake";
    md.shortName = "frlake";
    md.units = "1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "fraction_of_ice";
    md.shortName = "frseaice";
    md.units = "1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "isotropic_variance_of_filtered_topography";
    md.shortName = "varflt";
    md.units = "m+2";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "surface_velocity_scale";
    md.shortName = "ustar";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "surface_bouyancy_scale";
    md.shortName = "bstar";
    md.units = "ms-2";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "planetary_boundary_layer_height";
    md.shortName = "zpbl";
    md.units = "m";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "surface_exchange_coefficient_for_momentum";
    md.shortName = "cm";
    md.units = "kgm-2s-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "surface_exchange_coefficient_for_heat";
    md.shortName = "ct";
    md.units = "kgm-2s-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "surface_exchange_coefficient_for_moisture";
    md.shortName = "cq";
    md.units = "kgm-2s-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "KCBL_before_moist";
    md.shortName = "kcbl";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "surface_temp_before_moist";
    md.shortName = "tsm";
    md.units = "K";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "lower_index_where_Kh_greater_than_2";
    md.shortName = "khl";
    md.units = "1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "upper_index_where_Kh_greater_than_2";
    md.shortName = "khu";
    md.units = "1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "moist_air_density";
    md.shortName = "airdens";
    md.units = "kgm-3";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "zorl";
    md.shortName = "zorl";
    md.units = "cm";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "air_temperature_at_2m";
    md.shortName = "t2m";
    md.units = "K";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "tropopause_pressure";
    md.shortName = "tropopause_pressure";
    md.units = "Pa";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "geopotential_height";
    md.shortName = "geopotential_height";
    md.units = "m";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "geopotential_height_levels";
    md.shortName = "geopotential_height_levels";
    md.units = "m";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "height_above_mean_sea_level";
    md.shortName = "height_above_mean_sea_level";
    md.units = "m";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "height_above_mean_sea_level_at_surface";
    md.shortName = "height_above_mean_sea_level_at_surface";
    md.units = "m";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "geopotential_height_at_surface";
    md.shortName = "geopotential_height_at_surface";
    md.units = "m";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "water_vapor_mixing_ratio_wrt_dry_air";
    md.shortName = "water_vapor_mixing_ratio_wrt_dry_air";
    md.units = "1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "integrated_layer_ozone_in_air";
    md.shortName = "integrated_layer_ozone_in_air";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_content_of_cloud_liquid_water_in_atmosphere_layer";
    md.shortName = "mass_content_of_cloud_liquid_water_in_atmosphere_layer";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_content_of_cloud_ice_in_atmosphere_layer";
    md.shortName = "mass_content_of_cloud_ice_in_atmosphere_layer";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_content_of_rain_in_atmosphere_layer";
    md.shortName = "mass_content_of_rain_in_atmosphere_layer";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_content_of_snow_in_atmosphere_layer";
    md.shortName = "mass_content_of_snow_in_atmosphere_layer";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_content_of_graupel_in_atmosphere_layer";
    md.shortName = "mass_content_of_graupel_in_atmosphere_layer";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_content_of_hail_in_atmosphere_layer";
    md.shortName = "mass_content_of_hail_in_atmosphere_layer";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "effective_radius_of_cloud_liquid_water_particle";
    md.shortName = "effective_radius_of_cloud_liquid_water_particle";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "effective_radius_of_cloud_ice_particle";
    md.shortName = "effective_radius_of_cloud_ice_particle";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "effective_radius_of_rain_particle";
    md.shortName = "effective_radius_of_rain_particle";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "effective_radius_of_snow_particle";
    md.shortName = "effective_radius_of_snow_particle";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "effective_radius_of_graupel_particle";
    md.shortName = "effective_radius_of_graupel_particle";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "effective_radius_of_hail_particle";
    md.shortName = "effective_radius_of_hail_particle";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "totalSnowDepth_background_error";
    md.shortName = "totalSnowDepth_background_error";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "water_area_fraction";
    md.shortName = "water_area_fraction";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "land_area_fraction";
    md.shortName = "land_area_fraction";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "ice_area_fraction";
    md.shortName = "ice_area_fraction";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "surface_snow_area_fraction";
    md.shortName = "surface_snow_area_fraction";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "skin_temperature_at_surface_where_sea";
    md.shortName = "skin_temperature_at_surface_where_sea";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "skin_temperature_at_surface_where_land";
    md.shortName = "skin_temperature_at_surface_where_land";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "skin_temperature_at_surface_where_ice";
    md.shortName = "skin_temperature_at_surface_where_ice";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "skin_temperature_at_surface_where_snow";
    md.shortName = "skin_temperature_at_surface_where_snow";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "surface_snow_thickness";
    md.shortName = "surface_snow_thickness";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "vegetation_area_fraction";
    md.shortName = "vegetation_area_fraction";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "wind_speed_at_surface";
    md.shortName = "wind_speed_at_surface";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "wind_from_direction_at_surface";
    md.shortName = "wind_from_direction_at_surface";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "direction";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "leaf_area_index";
    md.shortName = "leaf_area_index";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "volume_fraction_of_condensed_water_in_soil";
    md.shortName = "volume_fraction_of_condensed_water_in_soil";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "soil_temperature";
    md.shortName = "soil_temperature";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "land_type_index_NPOESS";
    md.shortName = "land_type_index_NPOESS";
    md.units = "none";
    md.kind = "integer";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "vegetation_type_index";
    md.shortName = "vegetation_type_index";
    md.units = "none";
    md.kind = "integer";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "soil_type";
    md.shortName = "soil_type";
    md.units = "none";
    md.kind = "integer";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "surface_roughness_length";
    md.shortName = "surface_roughness_length";
    md.units = "m";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "wind_reduction_factor_at_10m";
    md.shortName = "wind_reduction_factor_at_10m";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "observable_domain_mask";
    md.shortName = "observable_domain_mask";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "surface_emissivity";
    md.shortName = "surface_emissivity";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "optical_thickness_of_atmosphere_layer";
    md.shortName = "optical_thickness_of_atmosphere_layer";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "toa_outgoing_radiance_per_unit_wavenumber";
    md.shortName = "toa_outgoing_radiance_per_unit_wavenumber";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "brightness_temperature";
    md.shortName = "brightness_temperature";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "brightness_temperature_assuming_clear_sky";
    md.shortName = "brightness_temperature_assuming_clear_sky";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "transmittances_of_atmosphere_layer";
    md.shortName = "transmittances_of_atmosphere_layer";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "weightingfunction_of_atmosphere_layer";
    md.shortName = "weightingfunction_of_atmosphere_layer";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "pressure_level_at_peak_of_weightingfunction";
    md.shortName = "pressure_level_at_peak_of_weightingfunction";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "average_surface_temperature_within_field_of_view";
    md.shortName = "average_surface_temperature_within_field_of_view";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "equivalent_reflectivity_factor";
    md.shortName = "equivalent_reflectivity_factor";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "sea_ice_category_area_fraction";
    md.shortName = "sea_ice_category_area_fraction";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "sea_ice_category_thickness";
    md.shortName = "sea_ice_category_thickness";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "sea_surface_height_above_geoid";
    md.shortName = "sea_surface_height_above_geoid";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "sea_water_potential_temperature";
    md.shortName = "sea_water_potential_temperature";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "sea_water_conservative_temperature";
    md.shortName = "sea_water_conservative_temperature";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "sea_water_absolute_salinity";
    md.shortName = "sea_water_absolute_salinity";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "sea_water_practical_salinity";
    md.shortName = "sea_water_practical_salinity";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "sea_water_salinity";
    md.shortName = "sea_water_salinity";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "sea_water_cell_thickness";
    md.shortName = "sea_water_cell_thickness";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "latent_heat_vaporization";
    md.shortName = "latent_heat_vaporization";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "net_downwelling_shortwave_radiation";
    md.shortName = "net_downwelling_shortwave_radiation";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "upward_latent_heat_flux_in_air";
    md.shortName = "upward_latent_heat_flux_in_air";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "upward_sensible_heat_flux_in_air";
    md.shortName = "upward_sensible_heat_flux_in_air";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "net_downwelling_longwave_radiation";
    md.shortName = "net_downwelling_longwave_radiation";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "friction_velocity_over_water";
    md.shortName = "friction_velocity_over_water";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    // Aerosols
    md.longName = "mass_fraction_of_dust001_in_air";
    md.shortName = "du001";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_dust002_in_air";
    md.shortName = "du002";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_dust003_in_air";
    md.shortName = "du003";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_dust004_in_air";
    md.shortName = "du004";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_dust005_in_air";
    md.shortName = "du005";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_sea_salt001_in_air";
    md.shortName = "ss001";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_sea_salt002_in_air";
    md.shortName = "ss002";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_sea_salt003_in_air";
    md.shortName = "ss003";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_sea_salt004_in_air";
    md.shortName = "ss004";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_sea_salt005_in_air";
    md.shortName = "ss005";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_hydrophobic_black_carbon_in_air";
    md.shortName = "bcphobic";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_hydrophilic_black_carbon_in_air";
    md.shortName = "bcphilic";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_hydrophobic_organic_carbon_in_air";
    md.shortName = "ocphobic";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_hydrophilic_organic_carbon_in_air";
    md.shortName = "ocphilic";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_nitrate001_in_air";
    md.shortName = "no3an1";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_nitrate002_in_air";
    md.shortName = "no3an2";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_nitrate003_in_air";
    md.shortName = "no3an3";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_so2_in_air";
    md.shortName = "so2";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mass_fraction_of_sulfate_in_air";
    md.shortName = "so4";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "volume_extinction_in_air_due_to_aerosol_particles_lambda1";
    md.shortName = "ext1";
    md.units = "km-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "volume_extinction_in_air_due_to_aerosol_particles_lambda2";
    md.shortName = "ext2";
    md.units = "km-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "volume_extinction_in_air_due_to_aerosol_particles_lambda3";
    md.shortName = "ext3";
    md.units = "km-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "pm25at";
    md.shortName = "pm25at";
    md.units = "none";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "pm25ac";
    md.shortName = "pm25ac";
    md.units = "none";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "pm25co";
    md.shortName = "pm25co";
    md.units = "none";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    // Trace gases
    md.longName = "volume_mixing_ratio_of_no2";
    md.shortName = "vmr_no2";
    md.units = "mol mol-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "volume_mixing_ratio_of_no";
    md.shortName = "vmr_no";
    md.units = "mol mol-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "volume_mixing_ratio_of_o3";
    md.shortName = "vmr_o3";
    md.units = "mol mol-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "volume_mixing_ratio_of_oh";
    md.shortName = "vmr_oh";
    md.units = "mol mol-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "volume_mixing_ratio_of_co";
    md.shortName = "vmr_co";
    md.units = "mol mol-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "form";
    md.shortName = "form";
    md.units = "none";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mole_fraction_of_carbon_dioxide_in_air";
    md.shortName = "co2";
    md.units = "none";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "ech4";
    md.shortName = "ech4";
    md.units = "none";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "ozone_mass_mixing_ratio";
    md.shortName = "o3mr";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "mole_fraction_of_ozone_in_air";
    md.shortName = "o3ppmv";
    md.units = "mole_fraction_of_ozone_in_air";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "odd_oxygen_mixing_ratio";
    md.shortName = "ox";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "full";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "raw_orography";
    md.shortName = "orog_raw";
    md.units = "m";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "filtered_orography";
    md.shortName = "orog_filt";
    md.units = "m";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";

    // Cold start variables
    addFieldMetadata(fieldsmetadata, nlev, md);
    md.longName = "d_grid_west_face_tangential_wind_component_cold_start";
    md.shortName = "u_w_cold";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "eastwest";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "d_grid_west_face_normal_wind_component_cold_start";
    md.shortName = "v_w_cold";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "eastwest";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "d_grid_south_face_tangential_wind_component_cold_start";
    md.shortName = "u_s_cold";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "northsouth";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "d_grid_south_face_normal_wind_component_cold_start";
    md.shortName = "v_s_cold";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "northsouth";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "u_component_of_native_D_grid_wind_cold_start";
    md.shortName = "ud_cold";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "northsouth";
    md.levels = "half";
    md.space = "vector";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "v_component_of_native_D_grid_wind_cold_start";
    md.shortName = "vd_cold";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "eastwest";
    md.levels = "half";
    md.space = "vector";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "surface_pressure_cold_start";
    md.shortName = "ps_cold";
    md.units = "Pa";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "1";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "air_temperature_cold_start";
    md.shortName = "t_cold";
    md.units = "K";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "air_pressure_thickness_cold_start";
    md.shortName = "delp_cold";
    md.units = "Pa";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "upward_air_velocity_cold_start";
    md.shortName = "w_cold";
    md.units = "ms-1";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "layer_thickness_cold_start";
    md.shortName = "zh_cold";
    md.units = "m";
    md.kind = "double";
    md.tracer = "false";
    md.horizontalStaggerLocation = "center";
    md.levels = "halfplusone";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "specific_humidity_cold_start";
    md.shortName = "sphum_cold";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "cloud_liquid_water_cold_start";
    md.shortName = "liq_wat_cold";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "ozone_mass_mixing_ratio_cold_start";
    md.shortName = "o3mr_cold";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "ice_wat_cold_start";
    md.shortName = "ice_wat_cold";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "rainwat_cold_start";
    md.shortName = "rainwat_cold";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "snowwat_cold_start";
    md.shortName = "snowwat_cold";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);

    md.longName = "graupel_cold_start";
    md.shortName = "graupel_cold";
    md.units = "kgkg-1";
    md.kind = "double";
    md.tracer = "true";
    md.horizontalStaggerLocation = "center";
    md.levels = "half";
    md.space = "magnitude";
    addFieldMetadata(fieldsmetadata, nlev, md);
  }
}  // namespace fv3jedi

