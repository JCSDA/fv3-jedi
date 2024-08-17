/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"

#include "fv3jedi/LinearVariableChange/LinearVariableChange.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/ModelData/ModelData.h"
#include "fv3jedi/State/State.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

LinearVariableChange::LinearVariableChange(const Geometry & geom,
                                           const eckit::Configuration & config)
  : geom_(geom), linearVariableChange_(), fieldsMetadata_(geom.fieldsMetaData()),
    vader_()
{
  params_.deserialize(config);
  eckit::LocalConfiguration variableChangeConfig = params_.toConfiguration();
  ModelData modelData{geom};
  eckit::LocalConfiguration vaderConfig;
  vaderConfig.set(vader::configCookbookKey,
                  variableChangeConfig.getSubConfiguration("vader custom cookbook"));
  vaderConfig.set(vader::configModelVarsKey, modelData.modelData());

  // Create vader with fv3-jedi custom cookbook
  vader_.reset(new vader::Vader(params_.linearVariableChangeParameters.value().vader,
                                vaderConfig));
}

// -------------------------------------------------------------------------------------------------

LinearVariableChange::~LinearVariableChange() {}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarTraj(const State & xfg, const oops::Variables & vars_out) {
  oops::Log::trace() << "LinearVariableChange::changeVarTraj starting" << std::endl;

  // Make sure vars are longname
  const oops::Variables vars = fieldsMetadata_.getLongNameFromAnyName(vars_out);

  // Call Vader's changeVarTraj to populate its trajectory
  State vader_xfg(xfg);

  // Record start variables
  oops::Variables varsFilled = vader_xfg.variablesIncludingInterfaceFields();

  oops::Variables varsVader = vars;
  varsVader -= varsFilled;  // Pass only the needed variables

  // Call Vader. On entry, varsVader holds the vars requested from Vader; on exit,
  // it holds the vars NOT fullfilled by Vader, i.e., the vars still to be requested elsewhere.
  // vader_.changeVarTraj also returns the variables fulfilled by Vader.
  atlas::FieldSet xfgfs;
  vader_xfg.toFieldSet(xfgfs);
  varsVaderPopulates_ = vader_->changeVarTraj(xfgfs, varsVader);
  if (varsVaderPopulates_.size() > 0) {
    varsFilled += varsVaderPopulates_;
    vader_xfg.updateFields(varsFilled);
    vader_xfg.fromFieldSet(xfgfs);
  }

  // Create the model variable change
  linearVariableChange_.reset(LinearVariableChangeFactory::create(vader_xfg, vader_xfg, geom_,
             params_.linearVariableChangeParameters.value()));
  oops::Log::trace() << "LinearVariableChange::changeVarTraj done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarTL(Increment & dx, const oops::Variables & vars_out) const {
  oops::Log::trace() << "LinearVariableChange::changeVarTL starting" << std::endl;
  // Make sure vars are longname
  const oops::Variables vars = fieldsMetadata_.getLongNameFromAnyName(vars_out);

  // If all variables already in incoming state just remove the no longer needed fields
  if (vars <= dx.variablesIncludingInterfaceFields()) {
    dx.updateFields(vars);
    oops::Log::trace() << "LinearVariableChange::changeVarTL done (identity)" << std::endl;
    return;
  }

  // Call Vader. On entry, varsVaderWillPopulate holds the vars requested from Vader; on exit,
  // it should be empty, since we know which variables Vader will do from the changeVarTraj
  // call.
  atlas::FieldSet dxfs;
  dx.toFieldSet(dxfs);
  oops::Variables varsVaderWillPopulate = varsVaderPopulates_;
  if (varsVaderWillPopulate.size() > 0) {
    vader_->changeVarTL(dxfs, varsVaderWillPopulate);
    ASSERT(varsVaderWillPopulate.size() == 0);

    // Set intermediate state for the Increment containing original fields plus the ones
    // Vader has done
    oops::Variables varsVader = dx.variablesIncludingInterfaceFields();
    varsVader += varsVaderPopulates_;
    dx.updateFields(varsVader);
    dx.fromFieldSet(dxfs);
  }

  // The to/fromFieldSet above is for a var change, so we know it's just adding/removing fields,
  // and is not editing the values within a particular field. So, the interface-specific fields are
  // still up to date (unless the var changes are coded incorrectly...).
  dx.setInterfaceFieldsOutOfDate(false);

  // Create output state
  Increment dxout(dx.geometry(), vars, dx.time());

  // Call fv3 linear variable change TL
  linearVariableChange_->multiply(dx, dxout);

  // Allocate any extra fields and remove fields no longer needed
  dx.updateFields(vars);

  // Copy data from temporary state
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::changeVarTL done" << dx << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarInverseTL(Increment & dx,
                                              const oops::Variables & vars_out,
                                              const bool force_varchange) const {
  oops::Log::trace() << "LinearVariableChange::changeVarInverseTL starting" << std::endl;

  // Make sure vars are longname
  const oops::Variables vars = fieldsMetadata_.getLongNameFromAnyName(vars_out);

  // If all variables already in incoming state just remove the no longer needed fields
  if ((vars <= dx.variablesIncludingInterfaceFields()) && !force_varchange) {
    dx.updateFields(vars);
    oops::Log::trace() << "LinearVariableChange::changeVarInverseTL done (identity)" << std::endl;
    return;
  }

  // Create output state
  Increment dxout(dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->multiplyInverse(dx, dxout);

  // Allocate any extra fields and remove fields no longer needed
  dx.updateFields(vars);

  // Copy data from temporary state
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::changeVarInverseTL done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarAD(Increment & dx, const oops::Variables & vars_out,
                                       const bool force_varchange) const {
  oops::Log::trace() << "LinearVariableChange::changeVarAD starting" << std::endl;
  // Make sure vars are longname
  const oops::Variables vars_long = fieldsMetadata_.getLongNameFromAnyName(vars_out);

  // If LVC is Model2X, then the adjoint needs to know about interface/model-specific variables
  // that oops isn't able to put into vars_out. Currently, only Model2GeoVaLs starts with model
  // vars, but this can also appear as "default" in the yaml...
  oops::Variables vars = vars_long;
  if (params_.linearVariableChangeParameters.value().name.value().value() == "Model2GeoVaLs"
      || params_.linearVariableChangeParameters.value().name.value().value() == "default") {
    // If geovals have winds, then we know model must have winds too. if "vars_out" has no winds,
    // that must be because there were only D-grid winds and these were hidden from OOPS.
    if (dx.variables().has("ua") || dx.variables().has("eastward_wind")) {
      if (!(vars_long.has("eastward_wind") || vars_long.has("u_component_of_native_D_grid_wind"))) {
        vars = oops::Variables((std::vector<std::string>){"u_component_of_native_D_grid_wind",
                                                          "v_component_of_native_D_grid_wind"});
        vars += vars_long;
      }
    }
  }

  // If all variables already in incoming state just remove the no longer needed fields
  if ((vars <= dx.variablesIncludingInterfaceFields()) && !force_varchange) {
    dx.updateFields(vars);
    oops::Log::trace() << "LinearVariableChange::changeVarAD done (identity)" << std::endl;
    return;
  }

  // Create dxin as a copy of dx, minus the variables created by Vader (in the forward direction)
  // This way we ensure the model code will not be able to do the adjoint for these vars
  Increment dxin(dx, true);  // true => full copy
  oops::Variables varsVaderDidntPopulate = dx.variablesIncludingInterfaceFields();
  varsVaderDidntPopulate -= varsVaderPopulates_;
  dxin.updateFields(varsVaderDidntPopulate);

  dx.updateFields(varsVaderPopulates_);
  // Create empty output state
  Increment dxout(dx.geometry(), vars, dx.time());

  // Call model's adjoint variable change.
  linearVariableChange_->multiplyAD(dxin, dxout);

  // dxout needs to temporarily have the variables that Vader populated put into it before
  // being passed into vader_.changeVarAD, so Vader can do its adjoints.
  atlas::FieldSet dxout_fs;
  dxout.toFieldSet(dxout_fs);
  oops::Variables varsVaderWillAdjoint = varsVaderPopulates_;
  if (varsVaderWillAdjoint.size() > 0) {
    atlas::FieldSet dx_fs;
    dx.toFieldSet(dx_fs);
    for (const auto field : dx_fs) {
      dxout_fs.add(field);
    }

    vader_->changeVarAD(dxout_fs, varsVaderWillAdjoint);

    // After changeVarAD, vader should have removed everything from varsVaderWillAdjoint,
    // indicating it did all the adjoints we expected it to.
    ASSERT(varsVaderWillAdjoint.size() == 0);
  }

  // Copy dxout into dx for return
  dx.updateFields(vars);
  dx.fromFieldSet(dxout_fs);

  // The to/fromFieldSet above is for a var change, so we know it's just adding/removing fields,
  // and is not editing the values within a particular field. So, the interface-specific fields are
  // still up to date (unless the var changes are coded incorrectly...).
  dx.setInterfaceFieldsOutOfDate(false);

  oops::Log::trace() << "LinearVariableChange::changeVarAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarInverseAD(Increment & dx,
                                              const oops::Variables & vars_out) const {
  oops::Log::trace() << "LinearVariableChange::changeVarInverseAD starting" << std::endl;

  // Make sure vars are longname
  const oops::Variables vars = fieldsMetadata_.getLongNameFromAnyName(vars_out);

  // If all variables already in incoming state just remove the no longer needed fields
  if (vars <= dx.variablesIncludingInterfaceFields()) {
    dx.updateFields(vars);
    oops::Log::trace() << "LinearVariableChange::changeVarInverseAD done (identity)" << std::endl;
    return;
  }

  // Create output state
  Increment dxout(dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->multiplyInverseAD(dx, dxout);

  // Allocate any extra fields and remove fields no longer needed
  dx.updateFields(vars);

  // Copy data from temporary state
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::changeVarInverseAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::print(std::ostream & os) const {
  os << "FV3-JEDI variable change";
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
