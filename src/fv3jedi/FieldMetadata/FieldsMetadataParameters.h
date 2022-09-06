/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#pragma once

#include <algorithm>
#include <iterator>
#include <map>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace fv3jedi {

// -----------------------------------------------------------------------------------------------

class FieldOverrideParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FieldOverrideParameters, Parameters)
 public:
  oops::RequiredParameter<std::string> longName{"long name", this};
  oops::OptionalParameter<std::string> varUnits{"units", this};
  oops::OptionalParameter<std::string> InOuName{"io name", this};
  oops::OptionalParameter<std::string> InOuFile{"io file", this};
  oops::OptionalParameter<std::string> IntrpTyp{"interpolation type", this};
  oops::OptionalParameter<std::string> IntrpMsk{"interpolation source point mask", this};
};

// -----------------------------------------------------------------------------------------------

class FieldsOverrideParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FieldsOverrideParameters, Parameters)
 public:
  oops::RequiredParameter<std::vector<FieldOverrideParameters>> fields{"field metadata", this};
};

// -----------------------------------------------------------------------------------------------

class FieldsMetadataParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FieldsMetadataParameters, Parameters)
 public:
  oops::OptionalParameter<std::string> override{"field metadata override", this};
};

// -----------------------------------------------------------------------------------------------

}  // namespace fv3jedi
