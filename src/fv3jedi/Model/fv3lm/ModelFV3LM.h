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

#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Model/fv3lm/ModelFV3LM.interface.h"
#include "fv3jedi/Model/ModelBase.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class ModelBias;
  class Increment;
  class State;

// -------------------------------------------------------------------------------------------------

class ModelFV3LM: public ModelBase,
                  private util::ObjectCounter<ModelFV3LM> {
 public:
  static const std::string classname() {return "fv3jedi::ModelFV3LM";}

  ModelFV3LM(const Geometry &, const eckit::Configuration &);
  ~ModelFV3LM();

/// Prepare model integration
  void initialize(State &) const override;

/// Model integration
  void step(State &, const ModelBias &) const override;

/// Finish model integration
  void finalize(State &) const override;

/// Utilities
  const util::Duration & timeResolution() const override {return tstep_;}

 private:
  void print(std::ostream &) const override;
  F90model keyConfig_;
  util::Duration tstep_;
  const Geometry geom_;
};
// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
