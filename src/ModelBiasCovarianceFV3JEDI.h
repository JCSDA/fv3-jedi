/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_MODELBIASCOVARIANCEFV3JEDI_H_
#define SRC_MODELBIASCOVARIANCEFV3JEDI_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace fv3jedi {
  class ModelBiasFV3JEDI;
  class ModelBiasIncrementFV3JEDI;
  class GeometryFV3JEDI;

// -----------------------------------------------------------------------------

class ModelBiasCovarianceFV3JEDI : public util::Printable,
                       private boost::noncopyable,
                       private util::ObjectCounter<ModelBiasCovarianceFV3JEDI> {
 public:
  static const std::string classname()
                                 {return "fv3jedi::ModelBiasCovarianceFV3JEDI";}

/// Constructor, destructor
  ModelBiasCovarianceFV3JEDI(const eckit::Configuration & conf,
                             const GeometryFV3JEDI &): conf_(conf) {}
  ~ModelBiasCovarianceFV3JEDI() {}

/// Linear algebra operators
  void linearize(const ModelBiasFV3JEDI &, const GeometryFV3JEDI &) {}
  void multiply(const ModelBiasIncrementFV3JEDI &,
                ModelBiasIncrementFV3JEDI) const {}
  void inverseMultiply(const ModelBiasIncrementFV3JEDI &,
                ModelBiasIncrementFV3JEDI) const {}
  void randomize(ModelBiasIncrementFV3JEDI &) const {}

  const eckit::Configuration & config() const {return conf_;}

 private:
  void print(std::ostream & os) const {}
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // SRC_MODELBIASCOVARIANCEFV3JEDI_H_
