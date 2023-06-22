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
