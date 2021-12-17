/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "boost/none_t.hpp"

#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/VariableChange/VariableChange.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

VariableChange::VariableChange(const Parameters_ & params, const Geometry & geometry)
  : fieldsMetadata_(geometry.fieldsMetaData()) {
  // Create the variable change
  variableChange_.reset(VariableChangeFactory::create(geometry,
                                                      params.variableChangeParameters.value()));
}

// -------------------------------------------------------------------------------------------------

VariableChange::~VariableChange() {}

// -------------------------------------------------------------------------------------------------

void VariableChange::changeVar(State & x, const oops::Variables & vars) const {
  // Trace
  oops::Log::trace() << "VariableChange::changeVar starting" << std::endl;

  // Convert vars to long names and field names
  const oops::Variables varsLongName = fieldsMetadata_.LongNameFromIONameLongNameOrFieldName(vars);

  // If all variables already in incoming state just remove the no longer needed fields
  if (varsLongName <= x.variablesLongName()) {
    x.updateFields(vars);
    oops::Log::info() << "VariableChange::changeVar done (identity)" << std::endl;
    return;
  }

  // Create output state
  State xout(*x.geometry(), vars, x.time());

  // Call variable change
  variableChange_->changeVar(x, xout);

  // Allocate any extra fields and remove fields no longer needed
  x.updateFields(vars);

  // Copy data from temporary state
  x = xout;

  // Trace
  oops::Log::trace() << "VariableChange::changeVar done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void VariableChange::changeVarInverse(State & x, const oops::Variables & vars) const {
  // Trace
  oops::Log::trace() << "VariableChange::changeVarInverse starting" << std::endl;

  // Convert vars to long names and field names
  const oops::Variables varsLongName = fieldsMetadata_.LongNameFromIONameLongNameOrFieldName(vars);

  // If all variables already in incoming state just remove the no longer needed fields
  if (varsLongName <= x.variablesLongName()) {
    x.updateFields(vars);
    oops::Log::info() << "VariableChange::changeVarInverse done (identity)" << std::endl;
    return;
  }

  // Create output state
  State xout(*x.geometry(), vars, x.time());

  // Call variable change
  variableChange_->changeVarInverse(x, xout);

  // Allocate any extra fields and remove fields no longer needed
  x.updateFields(vars);

  // Copy data from temporary state
  x = xout;

  // Trace
  oops::Log::trace() << "VariableChange::changeVarInverse done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void VariableChange::print(std::ostream & os) const {
  os << *variableChange_;
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
