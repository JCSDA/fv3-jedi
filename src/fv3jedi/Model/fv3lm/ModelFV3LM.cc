/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Model/fv3lm/ModelFV3LM.h"
#include "fv3jedi/ModelBias/ModelBias.h"
#include "fv3jedi/State/State.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------

/// Options taken by ModelFV3LM
class ModelFV3LMParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelFV3LMParameters, Parameters)

 public:
  oops::RequiredParameter<oops::Variables> modelVariables{ "model variables", this};
  oops::RequiredParameter<util::Duration> tstep{ "tstep", this};

  oops::RequiredParameter<int> lm_do_dyn{ "lm_do_dyn", this};
  oops::RequiredParameter<int> lm_do_trb{ "lm_do_trb", this};
  oops::RequiredParameter<int> lm_do_mst{ "lm_do_mst", this};

  oops::Parameter<bool> useInternalNamelist{ "use internal namelist", false, this};
  oops::OptionalParameter<std::string> namelistFilename{ "namelist filename", this};
};

// -------------------------------------------------------------------------------------------------
static oops::interface::ModelMaker<Traits, ModelFV3LM> makermodel_("FV3LM");
// -------------------------------------------------------------------------------------------------
ModelFV3LM::ModelFV3LM(const Geometry & resol, const eckit::Configuration & config)
  : keyConfig_(0), tstep_(0), geom_(resol), vars_(), an2model_(), finalVars_()
{
  oops::Log::trace() << "ModelFV3LM::ModelFV3LM starting" << std::endl;
  ModelFV3LMParameters params;
  params.deserialize(config);
  vars_ = oops::Variables(geom_.fieldsMetaData().getLongNameFromAnyName(params.modelVariables));
  tstep_ = util::Duration(config.getString("tstep"));
  fv3jedi_fv3lm_create_f90(config, geom_.toFortran(), keyConfig_);

  // This code prevents having to put VarChange params in yaml, since there is no actual user option
  eckit::LocalConfiguration varChangeConfig;
  varChangeConfig.set("variable change name", "Analysis2Model");
  an2model_.reset(new VariableChange(varChangeConfig, resol));

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
  if (!(vars_ <= xx.variablesIncludingInterfaceFields())) {
    finalVars_.reset(new oops::Variables(xx.variablesIncludingInterfaceFields()));
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
  fv3jedi_fv3lm_finalize_f90(keyConfig_, xx.toFortran());
  if (finalVars_) {
    const bool force_varchange = true;
    an2model_->changeVarInverse(xx, *finalVars_, force_varchange);
    finalVars_.reset(nullptr);  // reset to null for next initialize
  }
  oops::Log::trace() << "ModelFV3LM::finalize done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelFV3LM::print(std::ostream & os) const {
  os << "ModelFV3LM::print not implemented";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
