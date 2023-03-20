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

#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Utilities/Constants.h"
#include "fv3jedi/VariableChange/VariableChange.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

VariableChange::VariableChange(const Parameters_ & params, const Geometry & geometry)
  : fieldsMetadata_(geometry.fieldsMetaData()), vader_(), run_vader_(params.run_vader.value()),
    run_fv3jedi_(params.run_fv3jedi.value()) {
  // Create vader cookbook
  vader::cookbookConfigType vaderCustomCookbook =
                                        params.variableChangeParameters.value().vaderCustomCookbook;

  // Create the VaderConstructConfig object
  vader::VaderConstructConfig vaderConstructConfig(vaderCustomCookbook);

  // Add all constants to vader constructor config
  std::vector<std::string> allConstantsNames = getAllConstantsNames();
  for (std::string allConstantsName : allConstantsNames) {
    vaderConstructConfig.addToConfig<double>(allConstantsName, getConstant(allConstantsName));
  }

  // Add geometry data to vader constructor config
  vaderConstructConfig.addToConfig<double>("air_pressure_at_top_of_atmosphere_model",
                                           geometry.pTop());
  vaderConstructConfig.addToConfig<std::vector<double>>
                                  ("sigma_pressure_hybrid_coordinate_a_coefficient", geometry.ak());
  vaderConstructConfig.addToConfig<std::vector<double>>
                                  ("sigma_pressure_hybrid_coordinate_b_coefficient", geometry.bk());
  vaderConstructConfig.addToConfig<int>("nLevels", geometry.nLevels());

  // Create vader with fv3-jedi custom cookbook
  vader_.reset(new vader::Vader(params.variableChangeParameters.value().vader,
                                vaderConstructConfig));
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

  // Return if required vars in input
  // --------------------------------
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
  // vader_->changeVar also returns the variables fulfilled by Vader. These variables are allocated
  // and populated and added to the FieldSet (xfs).
  atlas::FieldSet xfs;
  x.toFieldSet(xfs);
  if (run_vader_) {
    varsFilled += vader_->changeVar(xfs, varsVader);
    x.updateFields(varsFilled);
  }
  x.fromFieldSet(xfs);

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

void VariableChange::changeVarInverse(State & x, const oops::Variables & vars_out) const {
  // Trace
  oops::Log::trace() << "VariableChange::changeVarInverse starting" << std::endl;

  // Make sure vars are longname
  // ---------------------------
  const oops::Variables vars = fieldsMetadata_.getLongNameFromAnyName(vars_out);

  // Return if required vars in input
  // --------------------------------
  if (vars <= x.variables()) {
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
  // vader_->changeVar also returns the variables fulfilled by Vader. These variables are allocated
  // and populated and added to the FieldSet (xfs).
  atlas::FieldSet xfs;
  x.toFieldSet(xfs);
  varsFilled += vader_->changeVar(xfs, varsVader);
  x.updateFields(varsFilled);
  x.fromFieldSet(xfs);

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
