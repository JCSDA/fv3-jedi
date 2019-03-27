/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_STATE_STATEFV3JEDI_H_
#define SRC_STATE_STATEFV3JEDI_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "GeometryFV3JEDI.h"
#include "IncrementFV3JEDI.h"
#include "StateFV3JEDIFortran.h"

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
  class GeometryFV3JEDI;
  class IncrementFV3JEDI;
  class GetValuesTrajFV3JEDI;

// FV3JEDI model state

// -----------------------------------------------------------------------------
class StateFV3JEDI : public util::Printable,
                private util::ObjectCounter<StateFV3JEDI> {
 public:
  static const std::string classname() {return "fv3jedi::StateFV3JEDI";}

// Constructor, destructor and basic operators
  StateFV3JEDI(const GeometryFV3JEDI &, const oops::Variables &,
               const util::DateTime &);
  StateFV3JEDI(const GeometryFV3JEDI &, const oops::Variables &,
               const eckit::Configuration &);
  StateFV3JEDI(const GeometryFV3JEDI &, const StateFV3JEDI &);
  StateFV3JEDI(const StateFV3JEDI &);
  virtual ~StateFV3JEDI();

  StateFV3JEDI & operator=(const StateFV3JEDI &);
  void zero();
  void accumul(const double &, const StateFV3JEDI &);

// Get state values at observation locations
  void getValues(const ufo::Locations &, const oops::Variables &,
                  ufo::GeoVaLs &) const;
  void getValues(const ufo::Locations &, const oops::Variables &,
                  ufo::GeoVaLs &, const GetValuesTrajFV3JEDI &) const;

// Interpolate state
  void changeResolution(const StateFV3JEDI & xx);

// Interactions with Increment
  StateFV3JEDI & operator+=(const IncrementFV3JEDI &);

// IO and diagnostics
  void read(const eckit::Configuration &);
  void analytic_init(const eckit::Configuration &, const GeometryFV3JEDI &);
  void write(const eckit::Configuration &) const;
  double norm() const;

// Utilities
  boost::shared_ptr<const GeometryFV3JEDI> geometry() const {return geom_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}

  int & toFortran() {return keyState_;}
  const int & toFortran() const {return keyState_;}

// Private methods and variables
 private:
  void print(std::ostream &) const;
  F90state keyState_;
  boost::shared_ptr<const GeometryFV3JEDI> geom_;
  oops::Variables vars_;
  util::DateTime time_;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // SRC_STATE_STATEFV3JEDI_H_
