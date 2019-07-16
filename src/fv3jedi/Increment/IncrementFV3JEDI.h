/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_INCREMENT_INCREMENTFV3JEDI_H_
#define FV3JEDI_INCREMENT_INCREMENTFV3JEDI_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "fv3jedi/Geometry/GeometryFV3JEDI.h"
#include "fv3jedi/Increment/IncrementFV3JEDIFortran.h"
#include "fv3jedi/State/StateFV3JEDI.h"
#include "oops/base/GeneralizedDepartures.h"
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
  class UnstructuredGrid;
}

namespace fv3jedi {
  class ModelBiasIncrementFV3JEDI;
  class ErrorCovarianceFV3JEDI;
  class StateFV3JEDI;
  class GetValuesTrajFV3JEDI;

// FV3JEDI increment

// -----------------------------------------------------------------------------

class IncrementFV3JEDI : public oops::GeneralizedDepartures,
                        public util::Printable,
                        private util::ObjectCounter<IncrementFV3JEDI> {
 public:
  static const std::string classname() {return "fv3jedi::IncrementFV3JEDI";}

/// Constructor, destructor
  IncrementFV3JEDI(const GeometryFV3JEDI &, const oops::Variables &,
                   const util::DateTime &);
  IncrementFV3JEDI(const GeometryFV3JEDI &, const IncrementFV3JEDI &);
  IncrementFV3JEDI(const IncrementFV3JEDI &, const bool);
  IncrementFV3JEDI(const IncrementFV3JEDI &);
  virtual ~IncrementFV3JEDI();

/// Basic operators
  void diff(const StateFV3JEDI &, const StateFV3JEDI &);
  void zero();
  void zero(const util::DateTime &);
  IncrementFV3JEDI & operator =(const IncrementFV3JEDI &);
  IncrementFV3JEDI & operator+=(const IncrementFV3JEDI &);
  IncrementFV3JEDI & operator-=(const IncrementFV3JEDI &);
  IncrementFV3JEDI & operator*=(const double &);
  void axpy(const double &, const IncrementFV3JEDI &, const bool check = true);
  double dot_product_with(const IncrementFV3JEDI &) const;
  void schur_product_with(const IncrementFV3JEDI &);
  void random();
  void dirac(const eckit::Configuration &);

/// Get increment values at observation locations
  void getValuesTL(const ufo::Locations &, const oops::Variables &,
                   ufo::GeoVaLs &, const GetValuesTrajFV3JEDI &) const;
  void getValuesAD(const ufo::Locations &, const oops::Variables &,
                   const ufo::GeoVaLs &, const GetValuesTrajFV3JEDI &);

/// Unstructured grid
  void ug_coord(oops::UnstructuredGrid &) const;
  void field_to_ug(oops::UnstructuredGrid &, const int &) const;
  void field_from_ug(const oops::UnstructuredGrid &, const int &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

  void updateTime(const util::Duration & dt) {time_ += dt;}

/// Other
  void accumul(const double &, const StateFV3JEDI &);
  void jnormgrad(const StateFV3JEDI &, const eckit::Configuration &);

// Utilities
  boost::shared_ptr<const GeometryFV3JEDI> geometry() const {return geom_;}

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
  boost::shared_ptr<const GeometryFV3JEDI> geom_;
  oops::Variables vars_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // FV3JEDI_INCREMENT_INCREMENTFV3JEDI_H_
