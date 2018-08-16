/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "src/ModelFV3JEDI.h"

#include "oops/util/Logger.h"
#include "ModelBiasFV3JEDI.h"
#include "FieldsFV3JEDI.h"
#include "Fortran.h"
#include "GeometryFV3JEDI.h"
#include "StateFV3JEDI.h"
#include "eckit/config/Configuration.h"
#include "oops/util/DateTime.h"
#include "UtilitiesFV3JEDI.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
ModelFV3JEDI::ModelFV3JEDI(const GeometryFV3JEDI & resol,
                            const eckit::Configuration & model)
  : keyConfig_(0), tstep_(0), geom_(resol),
  vars_(std::vector<std::string>{"u", "v", "t", "delp", "q"})
{
  oops::Log::trace() << "ModelFV3JEDI::ModelFV3JEDI" << std::endl;
  tstep_ = util::Duration(model.getString("tstep"));
  const eckit::Configuration * configc = &model;
  stageFv3Files(model);
  fv3jedi_model_setup_f90(&configc, geom_.toFortran(), keyConfig_);
  removeFv3Files();
  oops::Log::trace() << "ModelFV3JEDI created" << std::endl;
}
// -----------------------------------------------------------------------------
ModelFV3JEDI::~ModelFV3JEDI() {
  fv3jedi_model_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelFV3JEDI::initialize(StateFV3JEDI & xx) const {
  fv3jedi_model_prepare_integration_f90(keyConfig_, xx.fields().toFortran());
  oops::Log::debug() << "ModelFV3JEDI::initialize" << xx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void ModelFV3JEDI::step(StateFV3JEDI & xx, const ModelBiasFV3JEDI &) const {
  oops::Log::debug() << "ModelFV3JEDI::step fields in"
                     << xx.fields() << std::endl;
  fv3jedi_model_propagate_f90(keyConfig_, xx.fields().toFortran());
  xx.validTime() += tstep_;
  oops::Log::debug() << "ModelFV3JEDI::step fields out"
                     << xx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void ModelFV3JEDI::finalize(StateFV3JEDI & xx) const {
  oops::Log::debug() << "ModelFV3JEDI::finalize" << xx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
int ModelFV3JEDI::saveTrajectory(StateFV3JEDI & xx,
                                 const ModelBiasFV3JEDI &) const {
  int ftraj = 0;
  oops::Log::debug() << "ModelFV3JEDI::saveTrajectory fields in"
                     << xx.fields() << std::endl;
  fv3jedi_model_prop_traj_f90(keyConfig_, xx.fields().toFortran(), ftraj);
  ASSERT(ftraj != 0);
  oops::Log::debug() << "ModelFV3JEDI::saveTrajectory fields out"
                     << xx.fields() << std::endl;
  return ftraj;
}
// -----------------------------------------------------------------------------
void ModelFV3JEDI::print(std::ostream & os) const {
  os << "ModelFV3JEDI::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
