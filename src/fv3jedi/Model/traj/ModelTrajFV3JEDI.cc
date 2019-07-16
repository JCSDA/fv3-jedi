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

#include "fv3jedi/Geometry/GeometryFV3JEDI.h"
#include "fv3jedi/Model/traj/ModelTrajFV3JEDI.h"
#include "fv3jedi/ModelBias/ModelBiasFV3JEDI.h"
#include "fv3jedi/State/StateFV3JEDI.h"
#include "fv3jedi/Utilities/UtilitiesFV3JEDI.h"
#include "ModelTrajFV3JEDIFortran.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
static oops::ModelMaker<FV3JEDITraits, ModelTrajFV3JEDI> makermodel_("TRAJ");
// -----------------------------------------------------------------------------
ModelTrajFV3JEDI::ModelTrajFV3JEDI(const GeometryFV3JEDI & resol,
                            const eckit::Configuration & mconf)
  : keyConfig_(0), tstep_(0), geom_(resol), vars_(mconf)
{
  oops::Log::trace() << "ModelTrajFV3JEDI created" << std::endl;
}
// -----------------------------------------------------------------------------
ModelTrajFV3JEDI::~ModelTrajFV3JEDI() {
  oops::Log::trace() << "ModelTrajFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelTrajFV3JEDI::initialize(StateFV3JEDI & xx) const {
  oops::Log::debug() << "ModelTrajFV3JEDI::initialize" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelTrajFV3JEDI::step(StateFV3JEDI & xx, const ModelBiasFV3JEDI &) const {
  xx.validTime() += tstep_;
  oops::Log::debug() << "ModelTrajFV3JEDI::step" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelTrajFV3JEDI::finalize(StateFV3JEDI & xx) const {
  oops::Log::debug() << "ModelTrajFV3JEDI::finalize" << std::endl;
}
// -----------------------------------------------------------------------------
int ModelTrajFV3JEDI::saveTrajectory(StateFV3JEDI & xx,
                                 const ModelBiasFV3JEDI &) const {
  int ftraj = 0;
  fv3jedi_traj_prop_f90(keyConfig_, xx.toFortran(), ftraj);
  ASSERT(ftraj != 0);
  return ftraj;
}
// -----------------------------------------------------------------------------
void ModelTrajFV3JEDI::print(std::ostream & os) const {
  os << "ModelTrajFV3JEDI::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
