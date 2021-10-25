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

#include "oops/base/Variables.h"
#include "oops/interface/ModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/Geometry/Geometry.h"
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

class ModelPseudo: public oops::interface::ModelBase<Traits>,
                   private util::ObjectCounter<ModelPseudo> {
 public:
  static const std::string classname() {return "fv3jedi::ModelPseudo";}

  ModelPseudo(const Geometry &, const eckit::Configuration &);
  ~ModelPseudo();

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
  util::Duration tstep_;
  oops::Variables vars_;
  bool runstagecheck_;
  mutable bool runstage_ = true;
  std::unique_ptr<IOBase> io_;
};
// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
