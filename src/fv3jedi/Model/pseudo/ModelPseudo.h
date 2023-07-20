/*
 * (C) Copyright 2019-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/generic/ModelBase.h"
#include "oops/interface/ModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "fv3jedi/IO/Utils/IOBase.h"
#include "fv3jedi/Utilities/Traits.h"

namespace fv3jedi {
  class Geometry;
  class ModelBias;
  class Increment;
  class State;

// -------------------------------------------------------------------------------------------------

class ModelPseudoParameters : public oops::ModelParametersBase {
  OOPS_CONCRETE_PARAMETERS(ModelPseudoParameters, ModelParametersBase)
 public:
  oops::Parameter<bool> runstagecheck{ "run stage check", "turn off subsequent forecasts "
                                       "in multiple forecast applications such as outer loop data "
                                       "assimilation", false, this};
  oops::RequiredParameter<util::Duration> tstep{ "tstep", this};
  // Include IO parameters
  IOParametersWrapper ioParametersWrapper{this};
};

// -------------------------------------------------------------------------------------------------

class ModelPseudo: public oops::interface::ModelBase<Traits>,
                   private util::ObjectCounter<ModelPseudo> {
 public:
  static const std::string classname() {return "fv3jedi::ModelPseudo";}

  typedef ModelPseudoParameters Parameters_;

  ModelPseudo(const Geometry &, const Parameters_ &);
  ~ModelPseudo();

/// Prepare model integration
  void initialize(State &) const;

/// Model integration
  void step(State &, const ModelBias &) const;

/// Finish model integration
  void finalize(State &) const;

/// Utilities
  const util::Duration & timeResolution() const {return tstep_;}

 private:
  void print(std::ostream &) const;
  util::Duration tstep_;
  bool runstagecheck_;
  mutable bool runstage_ = true;
  std::unique_ptr<IOBase> io_;
};
// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
