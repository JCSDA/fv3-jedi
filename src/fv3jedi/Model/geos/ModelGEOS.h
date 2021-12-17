/*
 * (C) Copyright 2017 UCAR
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
#include "fv3jedi/Model/geos/ModelGEOS.interface.h"
#include "fv3jedi/Utilities/Traits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class ModelBias;
  class Increment;
  class State;

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/// Options taken by ModelGEOS
  class ModelGEOSParameters : public oops::ModelParametersBase {
    OOPS_CONCRETE_PARAMETERS(ModelGEOSParameters, ModelParametersBase)

   public:
    oops::RequiredParameter<oops::Variables> modelVariables{ "model variables", this};
    oops::RequiredParameter<std::string> geosRunDirectory{ "geos_run_directory", this};
    oops::RequiredParameter<util::Duration> tstep{ "tstep", this};

    oops::OptionalParameter<bool> reforecast{ "reforecast", this};
    oops::OptionalParameter<bool> esmfLogging{ "ESMF_Logging", this};
  };

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
/// FV3JEDI model definition.
/*!
 *  FV3JEDI nonlinear model definition and configuration parameters.
 */

class ModelGEOS: public oops::interface::ModelBase<Traits>,
                 private util::ObjectCounter<ModelGEOS> {
 public:
  typedef ModelGEOSParameters Parameters_;
  static const std::string classname() {return "fv3jedi::ModelGEOS";}

  ModelGEOS(const Geometry &, const Parameters_ &);
  ~ModelGEOS();

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
  char jedidir_[10000];
  char geosscrdir_[10000];
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
