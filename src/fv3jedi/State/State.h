/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_STATE_STATE_H_
#define FV3JEDI_STATE_STATE_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.interface.h"

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
  class Geometry;
  class Increment;

// FV3JEDI model state

// -------------------------------------------------------------------------------------------------
class State : public util::Printable, private util::ObjectCounter<State> {
 public:
  static const std::string classname() {return "fv3jedi::State";}

// Constructor, destructor and basic operators
  State(const Geometry &, const oops::Variables &, const util::DateTime &);
  State(const Geometry &, const eckit::Configuration &);
  State(const Geometry &, const State &);
  State(const State &);
  virtual ~State();

  State & operator=(const State &);
  void zero();
  void accumul(const double &, const State &);

// Interpolate state
  void changeResolution(const State & xx);

// Interactions with Increment
  State & operator+=(const Increment &);

// IO and diagnostics
  void read(const eckit::Configuration &);
  void analytic_init(const eckit::Configuration &, const Geometry &);
  void write(const eckit::Configuration &) const;
  double norm() const;

/// Serialize and deserialize
  size_t serialSize() const;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, size_t &);

// Utilities
  std::shared_ptr<const Geometry> geometry() const {return geom_;}
  const oops::Variables & variables() const {return vars_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}
  void updateTime(const util::Duration & dt) {time_ += dt;}

  int & toFortran() {return keyState_;}
  const int & toFortran() const {return keyState_;}

// Private methods and variables
 private:
  void print(std::ostream &) const;
  F90state keyState_;
  std::shared_ptr<const Geometry> geom_;
  oops::Variables vars_;
  util::DateTime time_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // FV3JEDI_STATE_STATE_H_
