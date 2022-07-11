/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/generic/ModelBase.h"
#include "oops/interface/ModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Model/fv3lm/ModelFV3LM.interface.h"
#include "fv3jedi/Utilities/Traits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class ModelBias;
  class Increment;
  class State;

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
/// Options taken by ModelFV3LM
  class ModelFV3LMParameters : public oops::ModelParametersBase {
    OOPS_CONCRETE_PARAMETERS(ModelFV3LMParameters, ModelParametersBase)

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
// -------------------------------------------------------------------------------------------------

class ModelFV3LM: public oops::interface::ModelBase<Traits>,
                  private util::ObjectCounter<ModelFV3LM> {
 public:
  typedef ModelFV3LMParameters Parameters_;

  static const std::string classname() {return "fv3jedi::ModelFV3LM";}

  ModelFV3LM(const Geometry &, const Parameters_ &);
  ~ModelFV3LM();

/// Prepare model integration
  void initialize(State &) const;

/// Model integration
  void step(State &, const ModelBias &) const;

/// Finish model integration
  void finalize(State &) const;

/// Utilities
  const util::Duration & timeResolution() const {return tstep_;}
  const oops::Variables & variables() const {return vars_;}

 private:
  void print(std::ostream &) const;
  F90model keyConfig_;
  util::Duration tstep_;
  const Geometry geom_;
  const oops::Variables vars_;
};
// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
