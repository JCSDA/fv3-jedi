/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "oops/util/Logger.h"

#include "fv3jedi/LinearVariableChange/LinearVariableChange.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

LinearVariableChange::LinearVariableChange(const Geometry & geom, const Parameters_ & params)
  : geom_(geom), params_(params), linearVariableChange_(),
    fieldsMetadata_(geom.fieldsMetaData()),
    vader_(params.linearVariableChangeParameters.value().vader) {}

// -------------------------------------------------------------------------------------------------

LinearVariableChange::~LinearVariableChange() {}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarTraj(const State & xfg, const oops::Variables & vars_out) {
  oops::Log::trace() << "LinearVariableChange::changeVarTraj starting" << std::endl;

  // Make sure vars are longname
  // ---------------------------
  const oops::Variables vars = fieldsMetadata_.getLongNameFromAnyName(vars_out);

  // Call Vader's changeVarTraj to populate its trajectory
  // ------------------------------------------------------
  State vader_xfg(xfg);

  // Record start variables
  oops::Variables varsPopulated = vader_xfg.variables();

  // Set first guess to have all possible variables
  oops::Variables varsTotal = vader_xfg.variables();
  varsTotal += vars;
  vader_xfg.updateFields(varsTotal);

  oops::Variables neededVars = vars;
  neededVars -= varsPopulated;  // Pass only the needed variables

  // Call Vader. On entry, neededVars holds the vars requested from Vader; on exit,
  // it holds the vars NOT fullfilled by Vader, i.e., the vars still to be requested elsewhere.
  // vader_.changeVarTraj also returns the variables fulfilled by Vader.
  atlas::FieldSet xfgfs;
  vader_xfg.toFieldSet(xfgfs);
  varsVaderPopulates_ = vader_.changeVarTraj(xfgfs, neededVars);
  varsPopulated += varsVaderPopulates_;
  vader_xfg.fromFieldSet(xfgfs);

  // Ahead of creating fv3jedi linear variable transform, add vader computed fields to first guess
  vader_xfg.updateFields(varsPopulated);

  // Create the model variable change
  linearVariableChange_.reset(LinearVariableChangeFactory::create(vader_xfg, vader_xfg, geom_,
             params_.linearVariableChangeParameters.value()));
  oops::Log::trace() << "LinearVariableChange::changeVarTraj done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarTL(Increment & dx, const oops::Variables & vars_out) const {
  // If all variables already in incoming state just remove the no longer needed fields
  if (vars_out <= dx.variables()) {
    dx.updateFields(vars_out);
    oops::Log::trace() << "LinearVariableChange::changeVarTL done (identity)" << std::endl;
    return;
  }

  // Make sure vars are longname
  // ---------------------------
  const oops::Variables vars = fieldsMetadata_.getLongNameFromAnyName(vars_out);

  // Record start variables
  oops::Variables varsPopulated = dx.variables();

  // Set first guess to have all possible variables
  oops::Variables varsVader = dx.variables();
  oops::Variables varsVaderWillPopulate = varsVaderPopulates_;
  varsVader += varsVaderWillPopulate;
  dx.updateFields(varsVader);

  // Call Vader. On entry, varsVaderWillPopulate holds the vars requested from Vader; on exit,
  // it should be empty, since we know which variables Vader will do from the changeVarTraj call.
  atlas::FieldSet dxfs;
  dx.toFieldSet(dxfs);

  vader_.changeVarTL(dxfs, varsVaderWillPopulate);
  ASSERT(varsVaderWillPopulate.size() == 0);
  dx.fromFieldSet(dxfs);

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

void LinearVariableChange::changeVarInverseTL(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::changeVarInverseTL starting" << std::endl;

  // If all variables already in incoming state just remove the no longer needed fields
  if (vars <= dx.variables()) {
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

void LinearVariableChange::changeVarAD(Increment & dx, const oops::Variables & vars_out) const {
  // If all variables already in incoming state just remove the no longer needed fields
  if (vars_out <= dx.variables()) {
    dx.updateFields(vars_out);
    oops::Log::trace() << "LinearVariableChange::changeVarAD done (identity)" << std::endl;
    return;
  }

  // Make sure vars are longname
  // ---------------------------
  const oops::Variables vars = fieldsMetadata_.getLongNameFromAnyName(vars_out);

  // Create dxin as a copy of dx, minus the variables created by Vader (in the forward direction)
  // This way we ensure the model code will not be able to do the adjoint for these vars
  Increment dxin(dx.geometry(), dx.variables(), dx.time());
  dxin = dx;
  oops::Variables varsVaderDidntPopulate = dx.variables();
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
  atlas::FieldSet dx_fs;
  dx.toFieldSet(dx_fs);
  for (const auto field : dx_fs) {
    dxout_fs.add(field);
  }

  oops::Variables varsVaderWillAdjoint = varsVaderPopulates_;
  vader_.changeVarAD(dxout_fs, varsVaderWillAdjoint);

  // After changeVarAD, vader should have removed everything from varsVaderWillAdjoint,
  // indicating it did all the adjoints we expected it to.
  ASSERT(varsVaderWillAdjoint.size() == 0);

  // Copy dxout into dx for return
  dx.updateFields(vars);
  dx.fromFieldSet(dxout_fs);

  oops::Log::trace() << "LinearVariableChange::changeVarAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarInverseAD(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::changeVarInverseAD starting" << std::endl;

  // If all variables already in incoming state just remove the no longer needed fields
  if (vars <= dx.variables()) {
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
