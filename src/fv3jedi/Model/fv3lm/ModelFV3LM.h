/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_MODEL_FV3LM_MODELFV3LM_H_
#define FV3JEDI_MODEL_FV3LM_MODELFV3LM_H_

#include <ostream>
#include <string>

#include "oops/base/ModelBase.h"
#include "oops/base/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
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

class ModelFV3LM: public oops::ModelBase<Traits>, private util::ObjectCounter<ModelFV3LM> {
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
#endif  // FV3JEDI_MODEL_FV3LM_MODELFV3LM_H_
