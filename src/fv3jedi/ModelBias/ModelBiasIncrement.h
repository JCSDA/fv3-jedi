/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <iostream>
#include <vector>

#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class ModelBias;
  class ModelBiasCovariance;
  class Geometry;

// -----------------------------------------------------------------------------

class ModelBiasIncrement : public util::Printable {
 public:
/// Constructor, destructor
  ModelBiasIncrement(const Geometry &,
                            const eckit::Configuration &) {}
  ModelBiasIncrement(const ModelBiasIncrement &,
                            const bool) {}
  ModelBiasIncrement(const ModelBiasIncrement &,
                            const eckit::Configuration &) {}
  ~ModelBiasIncrement() {}

/// Linear algebra operators
  void diff(const ModelBias &, const ModelBias &) {}
  void zero() {}
  ModelBiasIncrement & operator=(const
                                  ModelBiasIncrement &) {return *this;}
  ModelBiasIncrement & operator+=(const
                                  ModelBiasIncrement &) {return *this;}
  ModelBiasIncrement & operator-=(const
                                  ModelBiasIncrement &) {return *this;}
  ModelBiasIncrement & operator*=(const double) {return *this;}
  void axpy(const double, const ModelBiasIncrement &) {}
  double dot_product_with(const ModelBiasIncrement &)
                          const {return 0.0;}

/// Serialize-Deserialize
  size_t serialSize() {return 0;}
  void serialize(std::vector<double> & vect) const {}
  void deserialize(const std::vector<double> &, size_t & index) {}

/// I/O and diagnostics
  void read(const eckit::Configuration &) {}
  void write(const eckit::Configuration &) const {}
  double norm() const {return 0.0;}

 private:
  explicit ModelBiasIncrement(const ModelBiasCovariance &);
  void print(std::ostream & os) const {}
};

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
