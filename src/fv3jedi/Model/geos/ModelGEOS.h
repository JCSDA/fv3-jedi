/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Model/geos/ModelGEOS.interface.h"
#include "fv3jedi/Model/ModelBase.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class ModelBias;
  class Increment;
  class State;

// -----------------------------------------------------------------------------
/// FV3JEDI model definition.
/*!
 *  FV3JEDI nonlinear model definition and configuration parameters.
 */

class ModelGEOS: public ModelBase,
                 private util::ObjectCounter<ModelGEOS> {
 public:
  static const std::string classname() {return "fv3jedi::ModelGEOS";}

  ModelGEOS(const Geometry &, const eckit::Configuration &);
  ~ModelGEOS();

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
  char jedidir_[10000];
  char geosscrdir_[10000];
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
