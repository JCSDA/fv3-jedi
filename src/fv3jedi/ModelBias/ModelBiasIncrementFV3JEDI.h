/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_MODELBIAS_MODELBIASINCREMENTFV3JEDI_H_
#define FV3JEDI_MODELBIAS_MODELBIASINCREMENTFV3JEDI_H_

#include <iostream>

#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class ModelBiasFV3JEDI;
  class ModelBiasCovarianceFV3JEDI;
  class GeometryFV3JEDI;

// -----------------------------------------------------------------------------

class ModelBiasIncrementFV3JEDI : public util::Printable {
 public:
/// Constructor, destructor
  ModelBiasIncrementFV3JEDI(const GeometryFV3JEDI &,
                            const eckit::Configuration &) {}
  ModelBiasIncrementFV3JEDI(const ModelBiasIncrementFV3JEDI &,
                            const bool) {}
  ModelBiasIncrementFV3JEDI(const ModelBiasIncrementFV3JEDI &,
                            const eckit::Configuration &) {}
  ~ModelBiasIncrementFV3JEDI() {}

/// Linear algebra operators
  void diff(const ModelBiasFV3JEDI &, const ModelBiasFV3JEDI &) {}
  void zero() {}
  ModelBiasIncrementFV3JEDI & operator=(const
                                  ModelBiasIncrementFV3JEDI &) {return *this;}
  ModelBiasIncrementFV3JEDI & operator+=(const
                                  ModelBiasIncrementFV3JEDI &) {return *this;}
  ModelBiasIncrementFV3JEDI & operator-=(const
                                  ModelBiasIncrementFV3JEDI &) {return *this;}
  ModelBiasIncrementFV3JEDI & operator*=(const double) {return *this;}
  void axpy(const double, const ModelBiasIncrementFV3JEDI &) {}
  double dot_product_with(const ModelBiasIncrementFV3JEDI &)
                          const {return 0.0;}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0;}

 private:
  explicit ModelBiasIncrementFV3JEDI(const ModelBiasCovarianceFV3JEDI &);
  void print(std::ostream & os) const {}
};

// -----------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // FV3JEDI_MODELBIAS_MODELBIASINCREMENTFV3JEDI_H_
