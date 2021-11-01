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

class FieldParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FieldParameters, Parameters)
 public:
  oops::RequiredParameter<std::string> fieldName{ "FieldName", this};
  oops::Parameter<std::vector<std::string>> fieldIONames{ "FieldIONames", {}, this};
  oops::Parameter<std::string> interpType{ "InterpType", "default", this};
  oops::Parameter<std::string> ioFile{ "IOFile", "default", this};
  oops::Parameter<std::string> precisionKind{ "Kind", "double", this};
  oops::Parameter<std::string> levels{ "Levels", "full", this};
  oops::OptionalParameter<std::string> longName{ "LongName", this};
  oops::Parameter<std::string> space{ "Space", "magnitude", this};
  oops::Parameter<std::string> staggerLoc{ "StaggerLoc", "center", this};
  oops::Parameter<bool> tracer{ "Tracer", false, this};
  oops::RequiredParameter<std::string> units{ "Units", this};
};

// -----------------------------------------------------------------------------------------------

class FieldsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FieldsParameters, Parameters)
 public:
  oops::RequiredParameter<std::string> geometry{ "Geometry", this};
  oops::RequiredParameter<std::vector<FieldParameters>> fields{ "Fields", this};
};

// -----------------------------------------------------------------------------------------------

class FieldsetsParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FieldsetsParameters, Parameters)
 public:
  oops::RequiredParameter<std::string> fieldSet{ "fieldset", this};
};

// -----------------------------------------------------------------------------------------------

class FieldsMetadataParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(FieldsMetadataParameters, Parameters)
 public:
  oops::RequiredParameter<std::vector<FieldsetsParameters>> fieldSets{ "fieldsets", this};
};

// -----------------------------------------------------------------------------------------------

}  // namespace fv3jedi
