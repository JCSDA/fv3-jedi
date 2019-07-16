/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_MODEL_GEOS_MODELGEOSFV3JEDI_H_
#define FV3JEDI_MODEL_GEOS_MODELGEOSFV3JEDI_H_

#include <ostream>
#include <string>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "oops/base/ModelBase.h"
#include "oops/base/Variables.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/Geometry/GeometryFV3JEDI.h"
#include "fv3jedi/Utilities/FV3JEDITraits.h"
#include "ModelGEOSFV3JEDIFortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class ModelBiasFV3JEDI;
  class IncrementFV3JEDI;
  class StateFV3JEDI;

// -----------------------------------------------------------------------------
/// FV3JEDI model definition.
/*!
 *  FV3JEDI nonlinear model definition and configuration parameters.
 */

class ModelGEOSFV3JEDI: public oops::ModelBase<FV3JEDITraits>,
                        private util::ObjectCounter<ModelGEOSFV3JEDI> {
 public:
  static const std::string classname() {return "fv3jedi::ModelGEOSFV3JEDI";}

  ModelGEOSFV3JEDI(const GeometryFV3JEDI &, const eckit::Configuration &);
  ~ModelGEOSFV3JEDI();

/// Prepare model integration
  void initialize(StateFV3JEDI &) const;

/// Model integration
  void step(StateFV3JEDI &, const ModelBiasFV3JEDI &) const;
  int saveTrajectory(StateFV3JEDI &, const ModelBiasFV3JEDI &) const;

/// Finish model integration
  void finalize(StateFV3JEDI &) const;

/// Utilities
  const util::Duration & timeResolution() const {return tstep_;}
  const oops::Variables & variables() const {return vars_;}

 private:
  void print(std::ostream &) const;
  F90model keyConfig_;
  util::Duration tstep_;
  const GeometryFV3JEDI geom_;
  const oops::Variables vars_;
  char jedidir_[10000];
  char geosscrdir_[10000];
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_MODEL_GEOS_MODELGEOSFV3JEDI_H_
