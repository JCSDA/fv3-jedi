/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <string>
#include <vector>

namespace fv3jedi {

/// This static method is used to define the default value for the "vader custom cookbook" Parameter
/// for the VariableChange and LinearVariableChange classes.
// -------------------------------------------------------------------------------------------------

  static std::map<std::string, std::vector<std::string>> vaderFV3CustomCookbook() {
    return {
      // pt: from t and pkz, from t and ps
      {"air_potential_temperature",    {"AirPotentialTemperature_B", "AirPotentialTemperature_A"}},
      // P: from delp, from ps (and ak/bk)
      {"air_pressure_levels",          {"AirPressureAtInterface_B", "AirPressureAtInterface_A"}},
      // p: from pe
      {"air_pressure",                 {"AirPressure_A"}},
      // ln(p) from pe
      {"ln_air_pressure_at_interface", {"LnAirPressureAtInterface_A"}},
      // p^kappa from pe and ln(p)
      {"air_pressure_to_kappa",        {"AirPressureToKappa_A"}},
      // delp: from p
      {"air_pressure_thickness",       {"AirPressureThickness_A"}},
      // ps: from delp
      {"air_pressure_at_surface",      {"SurfaceAirPressure_A"}},
      // tv: from t and q
      {"virtual_temperature",          {"AirVirtualTemperature_A"}}
    };
  }

}  // namespace fv3jedi
