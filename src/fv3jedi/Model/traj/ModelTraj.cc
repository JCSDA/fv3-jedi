/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Model/traj/ModelTraj.h"
#include "fv3jedi/ModelBias/ModelBias.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Utilities/Utilities.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static oops::ModelMaker<Traits, ModelTraj> makermodel_("TRAJ");
// -------------------------------------------------------------------------------------------------
ModelTraj::ModelTraj(const Geometry & resol, const eckit::Configuration & mconf)
  : keyConfig_(0), tstep_(0), geom_(resol), vars_(mconf, "model variables")
{
  oops::Log::trace() << "ModelTraj created" << std::endl;
}
// -------------------------------------------------------------------------------------------------
ModelTraj::~ModelTraj() {
  oops::Log::trace() << "ModelTraj destructed" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelTraj::initialize(State & xx) const {
  oops::Log::debug() << "ModelTraj::initialize" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelTraj::step(State & xx, const ModelBias &) const {
  xx.validTime() += tstep_;
  oops::Log::debug() << "ModelTraj::step" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelTraj::finalize(State & xx) const {
  oops::Log::debug() << "ModelTraj::finalize" << std::endl;
}
// -------------------------------------------------------------------------------------------------
int ModelTraj::saveTrajectory(State & xx,
                                 const ModelBias &) const {
  int ftraj = 0;
  fv3jedi_traj_prop_f90(keyConfig_, xx.toFortran(), ftraj);
  ASSERT(ftraj != 0);
  return ftraj;
}
// -------------------------------------------------------------------------------------------------
void ModelTraj::print(std::ostream & os) const {
  os << "ModelTraj::print not implemented";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
