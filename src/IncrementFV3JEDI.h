/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3_JEDI_SRC_INCREMENTFV3JEDI_H_
#define FV3_JEDI_SRC_INCREMENTFV3JEDI_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "FieldsFV3JEDI.h"
#include "GeometryFV3JEDI.h"
#include "oops/base/GeneralizedDepartures.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/dot_product.h"

namespace eckit {
  class Configuration;
}

namespace ufo {
  class GeoVaLs;
}

namespace ioda {
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

/// Increment Class: Difference between two states
/*!
 *  Some fields that are present in a State may not be present in
 *  an Increment. The Increment contains everything that is needed by
 *  the tangent-linear and adjoint models.
 */

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
  void getValuesTL(const ioda::Locations &, const oops::Variables &,
                   ufo::GeoVaLs &, const GetValuesTrajFV3JEDI &) const;
  void getValuesAD(const ioda::Locations &, const oops::Variables &,
                   const ufo::GeoVaLs &, const GetValuesTrajFV3JEDI &);

/// Unstructured grid
  void ug_coord(oops::UnstructuredGrid &) const;
  void field_to_ug(oops::UnstructuredGrid &) const;
  void field_from_ug(const oops::UnstructuredGrid &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  const util::DateTime & validTime() const {return fields_->time();}
  util::DateTime & validTime() {return fields_->time();}
  void updateTime(const util::Duration & dt) {fields_->time() += dt;}

/// Access to fields
  FieldsFV3JEDI & fields() {return *fields_;}
  const FieldsFV3JEDI & fields() const {return *fields_;}

  boost::shared_ptr<const GeometryFV3JEDI> geometry() const {
    return fields_->geometry();
  }

/// Other
  void accumul(const double &, const StateFV3JEDI &);

/// Data
 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<FieldsFV3JEDI> fields_;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // FV3_JEDI_SRC_INCREMENTFV3JEDI_H_
