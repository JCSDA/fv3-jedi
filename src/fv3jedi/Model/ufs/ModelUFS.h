/*
 * (C) Copyright 2020 NOAA
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/interface/ModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Utilities/Traits.h"
#include "ModelUFS.interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class ModelBias;
  class Increment;
  class State;

// -------------------------------------------------------------------------------------------------

class ModelUFS: public oops::interface::ModelBase<Traits>,
                private util::ObjectCounter<ModelUFS> {
 public:
  static const std::string classname() {return "fv3jedi::ModelUFS";}

  ModelUFS(const Geometry &, const eckit::Configuration &);
  ~ModelUFS();

  void initialize(State &) const;
  void step(State &, const ModelBias &) const;

/// Finish model integration
  void finalize(State &) const;
  int saveTrajectory(State &, const ModelBias &) const;

  const util::Duration & timeResolution() const {return tstep_;}
  // note that this now excludes native u and v grid vars.
  // to include those, use variablesIncludingInterfaceFields()
  const oops::Variables & variables() const {return vars_;}

 private:
  void print(std::ostream &) const;
  F90model keyConfig_;
  util::Duration tstep_;
  util::Duration fclength_;
  const Geometry geom_;
  const oops::Variables vars_;
  char ufsdir_[10000];
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
