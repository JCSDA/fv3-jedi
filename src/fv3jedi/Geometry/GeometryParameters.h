/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "fv3jedi/FieldMetadata/FieldsMetadataParameters.h"
#include "fv3jedi/IO/Utils/IOBase.h"
#include "fv3jedi/State/State.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

class FMSinitParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FMSinitParameters, Parameters)

 public:
  oops::Parameter<std::string> fieldTableFilename{ "field table filename", "field_table", this};
  oops::Parameter<std::string> namelistFilename{ "namelist filename", "input.nml", this};
  oops::Parameter<int> stackmax{ "stackmax", 4000000, this};
};

// -------------------------------------------------------------------------------------------------

class GeometryParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(GeometryParameters, Parameters)

 public:
  oops::OptionalParameter<std::string> akbk{ "akbk", this};
  oops::Parameter<bool> doSchmidt{ "do_schmidt", false, this};
  oops::OptionalParameter<FMSinitParameters> fmsInit{ "fms initialization", this};
  oops::Parameter<bool> hydrostatic{ "hydrostatic", true, this};
  oops::Parameter<std::string> interpMethod{ "interpolation method", "barycent", this};
  oops::Parameter<std::vector<int>> ioLayout{ "io_layout", {1, 1}, this};
  oops::Parameter<std::vector<int>> layout{ "layout", {1, 1}, this};
  oops::Parameter<bool> logp{ "logp", false, this};
  oops::OptionalParameter<std::string> namelistFilename{"namelist filename", this};
  oops::Parameter<bool> nested{ "nested", false, this};
  oops::Parameter<int> ntiles{ "ntiles", 6, this};
  oops::OptionalParameter<int> npx{ "npx", this};
  oops::OptionalParameter<int> npy{ "npy", this};
  oops::OptionalParameter<int> npz{ "npz", this};
  oops::Parameter<int> iterator_dimension{ "iterator dimension", 2, this};
  oops::Parameter<int> nwat{ "nwat", 1, this};
  oops::OptionalParameter<StateParameters> orography{ "orography", this};
  oops::Parameter<bool> regional{ "regional", false, this};
  oops::Parameter<double> stretchFac{ "stretch_fac", 0.0, this};
  oops::Parameter<double> targetLat{ "target_lat", 0.0, this};
  oops::Parameter<double> targetLon{ "target_lon", 0.0, this};
  oops::Parameter<bool> useInternalNamelist{ "use internal namelist", false, this};
  oops::Parameter<bool> writeGeom{ "write geom", false, this};

  // Include FieldsMetadataParameters
  FieldsMetadataParameters fieldsMetadataParameters{this};
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
