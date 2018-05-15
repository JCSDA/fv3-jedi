/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef FV3JEDI_MODEL_STATEFV3JEDI_H_
#define FV3JEDI_MODEL_STATEFV3JEDI_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "FieldsFV3JEDI.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

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
  class GeometryFV3JEDI;
  class IncrementFV3JEDI;

/// FV3JEDI model state
/*!
 * A State contains everything that is needed to propagate the state
 * forward in time.
 */

// -----------------------------------------------------------------------------
class StateFV3JEDI : public util::Printable,
                private util::ObjectCounter<StateFV3JEDI> {
 public:
  static const std::string classname() {return "fv3jedi::StateFV3JEDI";}

/// Constructor, destructor
  StateFV3JEDI(const GeometryFV3JEDI &, const oops::Variables &, const util::DateTime &);  // Is it used?
  StateFV3JEDI(const GeometryFV3JEDI &, const eckit::Configuration &);
  StateFV3JEDI(const GeometryFV3JEDI &, const StateFV3JEDI &);
  StateFV3JEDI(const StateFV3JEDI &);
  virtual ~StateFV3JEDI();
  StateFV3JEDI & operator=(const StateFV3JEDI &);

/// Interpolate to observation location
  void interpolate(const ioda::Locations &, const oops::Variables &, ufo::GeoVaLs &) const;

/// Interpolate full fields
  void changeResolution(const StateFV3JEDI & xx);

/// Interactions with Increment
  StateFV3JEDI & operator+=(const IncrementFV3JEDI &);

/// Convert to/from generic unstructured grid
  void convert_to(oops::UnstructuredGrid &) const;
  void convert_from(const oops::UnstructuredGrid &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void analytic_init(const eckit::Configuration &, const GeometryFV3JEDI &);
  void write(const eckit::Configuration &) const;
  double norm() const {return fields_->norm();}
  const util::DateTime & validTime() const {return fields_->time();}
  util::DateTime & validTime() {return fields_->time();}

/// Access to fields
  FieldsFV3JEDI & fields() {return *fields_;}
  const FieldsFV3JEDI & fields() const {return *fields_;}

  boost::shared_ptr<const GeometryFV3JEDI> geometry() const {
    return fields_->geometry();
  }

/// Other
  void zero();
  void accumul(const double &, const StateFV3JEDI &);

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<FieldsFV3JEDI> fields_;
  boost::scoped_ptr<FieldsFV3JEDI> stash_;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // FV3JEDI_MODEL_STATEFV3JEDI_H_
