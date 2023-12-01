/*
 * (C) Copyright 2021-2023 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/VariableChangeParametersBase.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/ModelData/ModelData.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/VariableChange/VariableChange.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

VariableChange::VariableChange(const eckit::Configuration & config, const Geometry & geometry)
  : fieldsMetadata_(geometry.fieldsMetaData()), vader_(), run_vader_(), run_fv3jedi_()
{
  VariableChangeParametersWrapper params;
  params.deserialize(config);
  run_vader_ = params.run_vader.value();
  run_fv3jedi_ = params.run_fv3jedi.value();
  eckit::LocalConfiguration variableChangeConfig = params.toConfiguration();
  ModelData modelData{geometry};
  eckit::LocalConfiguration vaderConfig;
  vaderConfig.set(vader::configCookbookKey,
                  variableChangeConfig.getSubConfiguration("vader custom cookbook"));
  vaderConfig.set(vader::configModelVarsKey, modelData.modelData());

  // Create vader with fv3-jedi custom cookbook
  vader_.reset(new vader::Vader(params.variableChangeParameters.value().vader,
                                vaderConfig));
  // Create the variable change
  variableChange_.reset(VariableChangeFactory::create(geometry,
                                                      params.variableChangeParameters.value()));
}

// -------------------------------------------------------------------------------------------------

VariableChange::~VariableChange() {}

// -------------------------------------------------------------------------------------------------

void VariableChange::changeVar(State & x, const oops::Variables & vars_out) const {
  // Trace
  oops::Log::trace() << "VariableChange::changeVar starting" << std::endl;

  // Make sure vars are longname
  // ---------------------------
  const oops::Variables vars = fieldsMetadata_.getLongNameFromAnyName(vars_out);

  // Return if output vars already in input
  // --------------------------------------
  if (vars <= x.variables()) {
    x.updateFields(vars);
    oops::Log::info() << "VariableChange::changeVar done (identity)" << std::endl;
    return;
  }

  // Call Vader to perform first set of variable transforms
  // ------------------------------------------------------

  // Record start variables
  oops::Variables varsFilled = x.variables();

  oops::Variables varsVader = vars;
  varsVader -= varsFilled;  // Pass only the needed variables

  // Call Vader. On entry, varsVader holds the vars requested from Vader; on exit,
  // it holds the vars NOT fulfilled by Vader, i.e., the vars still to be requested elsewhere.
  // vader_->changeVar also returns the variables fulfilled by Vader. These variables are
  // allocated and populated and added to the FieldSet (xfs).
  if (run_vader_) {
    atlas::FieldSet xfs;
    x.toFieldSet(xfs);
    const oops::Variables varsVaderPopulated = vader_->changeVar(xfs, varsVader);
    if (varsVaderPopulated.size() > 0) {
      varsFilled += varsVaderPopulated;
      x.updateFields(varsFilled);
      x.fromFieldSet(xfs);
    }
  }

  // Perform fv3jedi factory variable change
  // ---------------------------------------

  // Create output state
  State xout(x.geometry(), vars, x.time());

  // Call variable change
  if (run_fv3jedi_) {
    variableChange_->changeVar(x, xout);
  }

  // Remove fields not in output
  x.updateFields(vars);

  // Copy data from temporary state
  x = xout;

  // Trace
  oops::Log::trace() << "VariableChange::changeVar done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void VariableChange::changeVarInverse(State & x, const oops::Variables & vars_out,
                                      const bool force_varchange) const {
  // Trace
  oops::Log::trace() << "VariableChange::changeVarInverse starting" << std::endl;

  // Make sure vars are longname
  // ---------------------------
  const oops::Variables vars = fieldsMetadata_.getLongNameFromAnyName(vars_out);

  // Return if output vars already in input
  // --------------------------------------
  if ((vars <= x.variables()) && !force_varchange) {
    x.updateFields(vars);
    oops::Log::info() << "VariableChange::changeVarInverse done (identity)" << std::endl;
    return;
  }

  // Call Vader to perform first set of variable transforms
  // ------------------------------------------------------

  // Record start variables
  oops::Variables varsFilled = x.variables();

  oops::Variables varsVader = vars;
  varsVader -= varsFilled;  // Pass only the needed variables

  // Call Vader. On entry, varsVader holds the vars requested from Vader; on exit,
  // it holds the vars NOT fulfilled by Vader, i.e., the vars still to be requested elsewhere.
  // vader_->changeVar also returns the variables fulfilled by Vader. These variables are
  // allocated and populated and added to the FieldSet (xfs).
  atlas::FieldSet xfs;
  x.toFieldSet(xfs);
  const oops::Variables varsVaderPopulated = vader_->changeVar(xfs, varsVader);
  if (varsVaderPopulated.size() > 0) {
    varsFilled += varsVaderPopulated;
    x.updateFields(varsFilled);
    x.fromFieldSet(xfs);
  }

  // Perform fv3jedi factory variable change
  // ---------------------------------------

  // Create output state
  State xout(x.geometry(), vars, x.time());

  // Call variable change
  variableChange_->changeVarInverse(x, xout);

  // Remove fields not in output
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
