/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef FV3JEDI_MODEL_FIELDSFV3JEDI_H_
#define FV3JEDI_MODEL_FIELDSFV3JEDI_H_

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "GeometryFV3JEDI.h"
#include "oops/base/Variables.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class UnstructuredGrid;
}

namespace ufo {
  class Locations;
  class GeoVaLs;
}

namespace fv3jedi {
// -----------------------------------------------------------------------------
/// Class to represent a FieldSet for the FV3JEDI model
class FieldsFV3JEDI : public util::Printable,
                    private util::ObjectCounter<FieldsFV3JEDI> {
 public:
  static const std::string classname() {return "fv3jedi::FieldsFV3JEDI";}

// Constructors and basic operators
  FieldsFV3JEDI(const GeometryFV3JEDI &, const oops::Variables &, const util::DateTime &);
  FieldsFV3JEDI(const FieldsFV3JEDI &, const GeometryFV3JEDI &);
  FieldsFV3JEDI(const FieldsFV3JEDI &, const oops::Variables &);
  FieldsFV3JEDI(const FieldsFV3JEDI &, const bool);
  FieldsFV3JEDI(const FieldsFV3JEDI &);
  ~FieldsFV3JEDI();

  void zero();
  void zero(const util::DateTime &);
  FieldsFV3JEDI & operator=(const FieldsFV3JEDI &);
  FieldsFV3JEDI & operator+=(const FieldsFV3JEDI &);
  FieldsFV3JEDI & operator-=(const FieldsFV3JEDI &);
  FieldsFV3JEDI & operator*=(const double &);
  void axpy(const double &, const FieldsFV3JEDI &);
  double dot_product_with(const FieldsFV3JEDI &) const;
  void schur_product_with(const FieldsFV3JEDI &);
  void random();
  void dirac(const eckit::Configuration &);

// Interpolate to given location
  void interpolate(const ufo::Locations &, const oops::Variables &, ufo::GeoVaLs &) const;
  void interpolateTL(const ufo::Locations &, const oops::Variables &, ufo::GeoVaLs &) const;
  void interpolateAD(const ufo::Locations &, const oops::Variables &, const ufo::GeoVaLs &);

// Interpolate full fields
  void changeResolution(const FieldsFV3JEDI &);
  void add(const FieldsFV3JEDI &);
  void diff(const FieldsFV3JEDI &, const FieldsFV3JEDI &);

// Convert to/from generic unstructured grid
  void convert_to(oops::UnstructuredGrid &) const;
  void convert_from(const oops::UnstructuredGrid &);

// Utilities
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  boost::shared_ptr<const GeometryFV3JEDI> geometry() const {return geom_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}

  int & toFortran() {return keyFlds_;}
  const int & toFortran() const {return keyFlds_;}

 private:
  void print(std::ostream &) const;
  F90flds keyFlds_;
  boost::shared_ptr<const GeometryFV3JEDI> geom_;
  oops::Variables vars_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_MODEL_FIELDSFV3JEDI_H_
