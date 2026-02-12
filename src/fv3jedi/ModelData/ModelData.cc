/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <ostream>
#include <string>
#include <vector>

#include "atlas/util/Config.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/ModelData/ModelData.h"
#include "fv3jedi/Utilities/Constants.h"

// -------------------------------------------------------------------------------------------------

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

ModelData::ModelData(const Geometry & geometry) :
    ak_(geometry.ak()),
    bk_(geometry.bk()),
    nLevels_(geometry.npz()),
    pTop_(geometry.pTop()) {}

// -------------------------------------------------------------------------------------------------

ModelData::~ModelData() {}

// -------------------------------------------------------------------------------------------------

const oops::Variables ModelData::defaultVariables() {
    return oops::Variables(std::vector<std::string>(
        {"air_temperature", "air_pressure", "air_pressure_levels",
         "water_area_fraction", "land_area_fraction", "ice_area_fraction",
         "surface_snow_area_fraction", "skin_temperature_at_surface_where_land",
         "skin_temperature_at_surface_where_ice", "skin_temperature_at_surface_where_snow",
         "skin_temperature_at_surface_where_sea", "vegetation_area_fraction", "leaf_area_index",
         "volume_fraction_of_condensed_water_in_soil", "soil_temperature", "surface_snow_thickness",
         "vegetation_type_index", "soil_type", "water_vapor_mixing_ratio_wrt_dry_air",
         "geopotential_height", "height_above_mean_sea_level_at_surface", "virtual_temperature",
         "mass_content_of_rain_in_atmosphere_layer", "mass_content_of_snow_in_atmosphere_layer",
         "mass_content_of_graupel_in_atmosphere_layer",
         "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
         "effective_radius_of_rain_particle", "effective_radius_of_snow_particle",
         "effective_radius_of_graupel_particle", "cloud_ice_number_concentration",
         "rain_number_concentration", "rain_water", "snow_water", "graupel",
         "mass_content_of_cloud_liquid_water_in_atmosphere_column",
         "mass_content_of_cloud_ice_in_atmosphere_column",
         "saturation_water_vapor_mixing_ratio_wrt_moist_air", "moist_air_density",
         "air_potential_temperature", "geopotential_height_at_surface",
         // TODO(AS): this has to be a variable that's derived from ocean and atmosphere,
         // needs to be in a CoupledVariableChange. For now, use the one computed in fv3-jedi
         "average_surface_temperature_within_field_of_view",
         "mole_fraction_of_ozone_in_air", "mole_fraction_of_carbon_dioxide_in_air",
         "effective_radius_of_cloud_liquid_water_particle", "land_type_index_NPOESS",
         "mass_content_of_cloud_ice_in_atmosphere_layer", "effective_radius_of_cloud_ice_particle",
         "wind_speed_at_surface", "wind_from_direction_at_surface", "tropopause_pressure",
         "eastward_wind", "northward_wind", "air_pressure_at_surface",
         "air_pressure_thickness", "water_vapor_mixing_ratio_wrt_moist_air",
         "cloud_liquid_ice", "cloud_liquid_water", "ozone_mass_mixing_ratio"}));
}

// -------------------------------------------------------------------------------------------------

const eckit::LocalConfiguration ModelData::modelData() const {
  eckit::LocalConfiguration modelData;

  // Add all constants to modelData config
  std::vector<std::string> allConstantsNames = getAllConstantsNames();
  for (std::string allConstantsName : allConstantsNames) {
    modelData.set(allConstantsName, getConstant(allConstantsName));
  }

  modelData.set("air_pressure_at_top_of_atmosphere_model", pTop_);
  modelData.set("sigma_pressure_hybrid_coordinate_a_coefficient", ak_);
  modelData.set("sigma_pressure_hybrid_coordinate_b_coefficient", bk_);
  modelData.set("nLevels", nLevels_);

  return modelData;
}

// -------------------------------------------------------------------------------------------------

void ModelData::print(std::ostream & os) const {
  os << "fv3jedi::ModelData::modelData(): " << modelData();
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
