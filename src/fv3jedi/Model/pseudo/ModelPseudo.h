/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_MODEL_PSEUDO_MODELPSEUDO_H_
#define FV3JEDI_MODEL_PSEUDO_MODELPSEUDO_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/interface/ModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Model/pseudo/ModelPseudo.interface.h"
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
/// FV3JEDI model definition.
/*!
 *  FV3JEDI nonlinear model definition and configuration parameters.
 */

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
  F90model keyConfig_;
  util::Duration tstep_;
  const Geometry geom_;
  const oops::Variables vars_;
  int runstagecheck_ = 0;
  mutable bool runstage_ = true;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_MODEL_PSEUDO_MODELPSEUDO_H_
