/*
 * (C) Copyright 2017-2020 UCAR
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
#include "fv3jedi/Model/fv3lm/ModelFV3LM.h"
#include "fv3jedi/ModelBias/ModelBias.h"
#include "fv3jedi/State/State.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static oops::interface::ModelMaker<Traits, ModelFV3LM> makermodel_("FV3LM");
// -------------------------------------------------------------------------------------------------
ModelFV3LM::ModelFV3LM(const Geometry & resol, const Parameters_ & params)
  : keyConfig_(0), tstep_(0), geom_(resol), vars_(params.modelVariables)
{
  oops::Log::trace() << "ModelFV3LM::ModelFV3LM" << std::endl;
  tstep_ = params.tstep;
  fv3jedi_fv3lm_create_f90(params.toConfiguration(), geom_.toFortran(), keyConfig_);
  oops::Log::trace() << "ModelFV3LM created" << std::endl;
}
// -------------------------------------------------------------------------------------------------
ModelFV3LM::~ModelFV3LM() {
  fv3jedi_fv3lm_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelFV3LM destructed" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelFV3LM::initialize(State & xx) const {
  fv3jedi_fv3lm_initialize_f90(keyConfig_, xx.toFortran());
  oops::Log::debug() << "ModelFV3LM::initialize" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelFV3LM::step(State & xx, const ModelBias &) const {
  xx.validTime() += tstep_;
  fv3jedi_fv3lm_step_f90(keyConfig_, xx.toFortran(), geom_.toFortran());
  oops::Log::debug() << "ModelFV3LM::step" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelFV3LM::finalize(State & xx) const {
  fv3jedi_fv3lm_finalize_f90(keyConfig_, xx.toFortran());
  oops::Log::debug() << "ModelFV3LM::finalize" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelFV3LM::print(std::ostream & os) const {
  os << "ModelFV3LM::print not implemented";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
