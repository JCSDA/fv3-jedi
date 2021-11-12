/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <fstream>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/GetValues/LinearGetValues.interface.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"

// -------------------------------------------------------------------------------------------------

namespace ufo {
  class GeoVaLs;
  class Locations;
}

namespace fv3jedi {
  class State;
  class Geometry;

// -------------------------------------------------------------------------------------------------

class LinearGetValues : public util::Printable, private util::ObjectCounter<LinearGetValues> {
 public:
  static const std::string classname() {return "fv3jedi::LinearGetValues";}

  LinearGetValues(const Geometry &, const ufo::Locations &, const eckit::Configuration &);
  virtual ~LinearGetValues();

  void setTrajectory(const State & state, const util::DateTime & t1, const util::DateTime & t2,
                     ufo::GeoVaLs & geovals);
  void fillGeoVaLsTL(const Increment & inc, const util::DateTime & t1, const util::DateTime & t2,
                     ufo::GeoVaLs & geovals) const;
  void fillGeoVaLsAD(Increment & inc, const util::DateTime & t1, const util::DateTime & t2,
                     const ufo::GeoVaLs & geovals) const;

 private:
  void print(std::ostream &) const;
  F90lineargetvalues keyLinearGetValues_;
  ufo::Locations locs_;
  std::shared_ptr<const Geometry> geom_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
