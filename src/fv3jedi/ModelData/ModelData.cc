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
    nLevels_(geometry.nLevels()),
    pTop_(geometry.pTop()) {}

// -------------------------------------------------------------------------------------------------

ModelData::~ModelData() {}

// -------------------------------------------------------------------------------------------------

const oops::Variables ModelData::defaultVariables() {
    return oops::Variables(std::vector<std::string>(
        {"air_temperature", "air_pressure", "air_pressure_levels",
         "water_area_fraction", "land_area_fraction", "ice_area_fraction",
         "surface_snow_area_fraction", "surface_temperature_where_land",
         "surface_temperature_where_ice", "surface_temperature_where_snow",
         "surface_temperature_where_sea", "vegetation_area_fraction", "leaf_area_index",
         "volume_fraction_of_condensed_water_in_soil", "soil_temperature", "surface_snow_thickness",
         "vegetation_type_index", "soil_type", "humidity_mixing_ratio",
         "mole_fraction_of_ozone_in_air", "mole_fraction_of_carbon_dioxide_in_air",
         "mass_content_of_cloud_liquid_water_in_atmosphere_layer",
         "effective_radius_of_cloud_liquid_water_particle",
         "mass_content_of_cloud_ice_in_atmosphere_layer", "effective_radius_of_cloud_ice_particle",
         "surface_wind_speed", "surface_wind_from_direction"}));
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
