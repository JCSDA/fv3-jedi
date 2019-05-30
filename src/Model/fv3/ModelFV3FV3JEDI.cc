/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "ModelFV3FV3JEDIFortran.h"
#include "GeometryFV3JEDI.h"
#include "ModelBiasFV3JEDI.h"
#include "src/Model/fv3/ModelFV3FV3JEDI.h"
#include "StateFV3JEDI.h"
#include "UtilitiesFV3JEDI.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
static oops::ModelMaker<FV3JEDITraits, ModelFV3FV3JEDI> makermodel_("FV3");
// -----------------------------------------------------------------------------
ModelFV3FV3JEDI::ModelFV3FV3JEDI(const GeometryFV3JEDI & resol,
                            const eckit::Configuration & mconf)
  : keyConfig_(0), tstep_(0), geom_(resol), vars_(mconf)
{
  oops::Log::trace() << "ModelFV3FV3JEDI::ModelFV3FV3JEDI" << std::endl;
  tstep_ = util::Duration(mconf.getString("tstep"));
  const eckit::Configuration * configc = &mconf;
  stageFv3Files(mconf);
  fv3jedi_fv3_create_f90(&configc, geom_.toFortran(), keyConfig_);
  removeFv3Files();
  oops::Log::trace() << "ModelFV3FV3JEDI created" << std::endl;
}
// -----------------------------------------------------------------------------
ModelFV3FV3JEDI::~ModelFV3FV3JEDI() {
  fv3jedi_fv3_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelFV3FV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelFV3FV3JEDI::initialize(StateFV3JEDI & xx) const {
  fv3jedi_fv3_initialize_f90(keyConfig_, xx.toFortran());
  oops::Log::debug() << "ModelFV3FV3JEDI::initialize" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelFV3FV3JEDI::step(StateFV3JEDI & xx, const ModelBiasFV3JEDI &) const {
  xx.validTime() += tstep_;
  util::DateTime * dtp = &xx.validTime();
  fv3jedi_fv3_step_f90(keyConfig_, xx.toFortran(), geom_.toFortran(), &dtp);
  oops::Log::debug() << "ModelFV3FV3JEDI::step" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelFV3FV3JEDI::finalize(StateFV3JEDI & xx) const {
  fv3jedi_fv3_finalize_f90(keyConfig_, xx.toFortran());
  oops::Log::debug() << "ModelFV3FV3JEDI::finalize" << std::endl;
}
// -----------------------------------------------------------------------------
int ModelFV3FV3JEDI::saveTrajectory(StateFV3JEDI & xx,
                                 const ModelBiasFV3JEDI &) const {
  ABORT("Model FV3 should not be used for the trajectory");
}
// -----------------------------------------------------------------------------
void ModelFV3FV3JEDI::print(std::ostream & os) const {
  os << "ModelFV3FV3JEDI::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
