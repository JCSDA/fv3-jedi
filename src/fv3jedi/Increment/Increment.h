/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_INCREMENT_INCREMENT_H_
#define FV3JEDI_INCREMENT_INCREMENT_H_

#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.interface.h"
#include "fv3jedi/State/State.h"
#include "oops/base/GeneralizedDepartures.h"
#include "oops/base/LocalIncrement.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace ufo {
  class GeoVaLs;
  class Locations;
}

namespace oops {
  class Variables;
}

namespace fv3jedi {
  class ModelBiasIncrement;
  class ErrorCovariance;
  class State;

// FV3JEDI increment

// -----------------------------------------------------------------------------

class Increment : public oops::GeneralizedDepartures,
                        public util::Printable,
                        private util::ObjectCounter<Increment> {
 public:
  static const std::string classname() {return "fv3jedi::Increment";}

/// Constructor, destructor
  Increment(const Geometry &, const oops::Variables &,
                   const util::DateTime &);
  Increment(const Geometry &, const Increment &);
  Increment(const Increment &, const bool);
  Increment(const Increment &);
  virtual ~Increment();

/// Basic operators
  void diff(const State &, const State &);
  void zero();
  void zero(const util::DateTime &);
  Increment & operator =(const Increment &);
  Increment & operator+=(const Increment &);
  Increment & operator-=(const Increment &);
  Increment & operator*=(const double &);
  void axpy(const double &, const Increment &, const bool check = true);
  double dot_product_with(const Increment &) const;
  void schur_product_with(const Increment &);
  void random();
  void dirac(const eckit::Configuration &);

/// Get/Set increment values at grid points
  oops::LocalIncrement getLocal(const GeometryIterator &) const;
  void setLocal(const oops::LocalIncrement &, const GeometryIterator &);

/// ATLAS
  void setAtlas(atlas::FieldSet *) const;
  void toAtlas(atlas::FieldSet *) const;
  void fromAtlas(atlas::FieldSet *);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

  void updateTime(const util::Duration & dt) {time_ += dt;}

/// Other
  void accumul(const double &, const State &);
  void jnormgrad(const State &, const eckit::Configuration &);

/// Serialize and deserialize
  size_t serialSize() const;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, size_t &);

// Utilities
  boost::shared_ptr<const Geometry> geometry() const {return geom_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}

  int & toFortran() {return keyInc_;}
  const int & toFortran() const {return keyInc_;}

// Private methods and variables
 private:
  void print(std::ostream &) const;
  F90inc keyInc_;
  boost::shared_ptr<const Geometry> geom_;
  oops::Variables vars_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // FV3JEDI_INCREMENT_INCREMENT_H_
