/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_MODELBIAS_MODELBIASFV3JEDI_H_
#define FV3JEDI_MODELBIAS_MODELBIASFV3JEDI_H_

#include <iostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class GeometryFV3JEDI;
  class ModelBiasIncrementFV3JEDI;

/// Model error for the FV3JEDI model.
/*!
 * This class is used to manipulate parameters of the model that
 * can be estimated in the assimilation. This includes model bias for
 * example but could be used for other parameters to be estimated.
 * This is sometimes referred to as augmented state or augmented
 * control variable in the litterature.
 * The augmented state is understood here as an augmented 4D state.
 */

// -----------------------------------------------------------------------------

class ModelBiasFV3JEDI : public util::Printable,
                       private boost::noncopyable,
                       private util::ObjectCounter<ModelBiasFV3JEDI> {
 public:
  static const std::string classname() {return "fv3jedi::ModelBiasFV3JEDI";}

  ModelBiasFV3JEDI(const GeometryFV3JEDI &, const eckit::Configuration &) {}
  ModelBiasFV3JEDI(const GeometryFV3JEDI &, const ModelBiasFV3JEDI &) {}
  ModelBiasFV3JEDI(const ModelBiasFV3JEDI &, const bool) {}
  ~ModelBiasFV3JEDI() {}

  ModelBiasFV3JEDI & operator+=(const
                          ModelBiasIncrementFV3JEDI &) {return *this;}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0;}

 private:
  void print(std::ostream & os) const {}
};

// -----------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // FV3JEDI_MODELBIAS_MODELBIASFV3JEDI_H_
