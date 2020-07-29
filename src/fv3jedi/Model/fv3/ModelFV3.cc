/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Model/fv3/ModelFV3.h"
#include "fv3jedi/ModelBias/ModelBias.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Utilities/Utilities.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static oops::ModelMaker<Traits, ModelFV3> makermodel_("FV3");
// -------------------------------------------------------------------------------------------------
ModelFV3::ModelFV3(const Geometry & resol, const eckit::Configuration & mconf)
  : keyConfig_(0), tstep_(0), geom_(resol), vars_(mconf, "model variables")
{
  oops::Log::trace() << "ModelFV3::ModelFV3" << std::endl;
  tstep_ = util::Duration(mconf.getString("tstep"));
  const eckit::Configuration * configc = &mconf;
  stageFv3Files(mconf, geom_.getComm());
  fv3jedi_fv3_create_f90(&configc, geom_.toFortran(), keyConfig_);
  removeFv3Files(geom_.getComm());
  oops::Log::trace() << "ModelFV3 created" << std::endl;
}
// -------------------------------------------------------------------------------------------------
ModelFV3::~ModelFV3() {
  fv3jedi_fv3_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelFV3 destructed" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelFV3::initialize(State & xx) const {
  util::DateTime * dtp = &xx.validTime();
  fv3jedi_fv3_initialize_f90(keyConfig_, xx.toFortran(), &dtp);
  oops::Log::debug() << "ModelFV3::initialize" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelFV3::step(State & xx, const ModelBias &) const {
  util::DateTime * dtp = &xx.validTime();
  fv3jedi_fv3_step_f90(keyConfig_, xx.toFortran(), &dtp);
  xx.validTime() += tstep_;
  oops::Log::debug() << "ModelFV3::step" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelFV3::finalize(State & xx) const {
  fv3jedi_fv3_finalize_f90(keyConfig_, xx.toFortran());
  oops::Log::debug() << "ModelFV3::finalize" << std::endl;
}
// -------------------------------------------------------------------------------------------------
int ModelFV3::saveTrajectory(State & xx,
                                 const ModelBias &) const {
  ABORT("Model FV3 should not be used for the trajectory");
  abort(); /* Prevent g++ missing return statement warning */
}
// -------------------------------------------------------------------------------------------------
void ModelFV3::print(std::ostream & os) const {
  os << "ModelFV3::print not implemented";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
