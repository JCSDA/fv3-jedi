/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/interface/ModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Model/fv3lm/ModelFV3LM.interface.h"
#include "fv3jedi/Utilities/Traits.h"
#include "fv3jedi/VariableChange/Base/VariableChangeBase.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class ModelBias;
  class Increment;
  class State;

// -------------------------------------------------------------------------------------------------

class ModelFV3LM: public oops::interface::ModelBase<Traits>,
                  private util::ObjectCounter<ModelFV3LM> {
 public:
  static const std::string classname() {return "fv3jedi::ModelFV3LM";}

  ModelFV3LM(const Geometry &, const eckit::Configuration &);
  ~ModelFV3LM();

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
  F90model keyConfig_;
  util::Duration tstep_;
  const Geometry geom_;
  oops::Variables vars_;
  std::unique_ptr<VariableChange> an2model_;
  mutable std::unique_ptr<const oops::Variables> finalVars_;
};
// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
