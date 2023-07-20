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
  : keyConfig_(0), tstep_(0), geom_(resol),
    vars_(geom_.fieldsMetaData().getLongNameFromAnyName(params.modelVariables)),
    an2model_(), finalVars_()
{
  oops::Log::trace() << "ModelFV3LM::ModelFV3LM starting" << std::endl;
  tstep_ = params.tstep;

  // This code prevents having to put VarChange params in yaml, since there is no actual user option
  VariableChangeParametersWrapper varChangeParams;
  eckit::LocalConfiguration varChangeConfig;
  varChangeConfig.set("variable change name", "Analysis2Model");
  varChangeParams.deserialize(varChangeConfig);
  an2model_.reset(new VariableChange(varChangeParams, resol));

  fv3jedi_fv3lm_create_f90(params.toConfiguration(), geom_.toFortran(), keyConfig_);
  oops::Log::trace() << "ModelFV3LM::ModelFV3LM done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
ModelFV3LM::~ModelFV3LM() {
  oops::Log::trace() << "ModelFV3LM::~ModelFV3LM starting" << std::endl;
  fv3jedi_fv3lm_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelFV3LM::~ModelFV3LM done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelFV3LM::initialize(State & xx) const {
  oops::Log::trace() << "ModelFV3LM::initialize starting" << std::endl;
  ASSERT_MSG(!finalVars_, "finalVars_ should always be null when calling initialize");
  if (!(vars_ <= xx.variables())) {
    finalVars_.reset(new oops::Variables(xx.variables()));
    an2model_->changeVar(xx, vars_);
  }
  fv3jedi_fv3lm_initialize_f90(keyConfig_, xx.toFortran());
  oops::Log::trace() << "ModelFV3LM::initialize done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelFV3LM::step(State & xx, const ModelBias &) const {
  oops::Log::trace() << "ModelFV3LM::step starting" << std::endl;
  xx.validTime() += tstep_;
  fv3jedi_fv3lm_step_f90(keyConfig_, xx.toFortran(), geom_.toFortran());
  oops::Log::trace() << "ModelFV3LM::step done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelFV3LM::finalize(State & xx) const {
  oops::Log::trace() << "ModelFV3LM::finalize starting" << std::endl;
  if (finalVars_) {
    an2model_->changeVarInverse(xx, *finalVars_);
    finalVars_.reset(nullptr);  // reset to null for next initialize
  }
  fv3jedi_fv3lm_finalize_f90(keyConfig_, xx.toFortran());
  oops::Log::trace() << "ModelFV3LM::finalize done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelFV3LM::print(std::ostream & os) const {
  os << "ModelFV3LM::print not implemented";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
