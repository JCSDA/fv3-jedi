/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

#include "fv3jedi/GetValues/GetValues.interface.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "ufo/Locations.h"

// -------------------------------------------------------------------------------------------------

namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace ufo {
  class GeoVaLs;
}

namespace fv3jedi {
  class State;
  class Geometry;

// -------------------------------------------------------------------------------------------------

class GetValues : public util::Printable, private util::ObjectCounter<GetValues> {
 public:
  static const std::string classname() {return "fv3jedi::GetValues";}

  GetValues(const Geometry &, const ufo::Locations &, const eckit::Configuration &);
  virtual ~GetValues();

  void fillGeoVaLs(const State &, const util::DateTime &, const util::DateTime &,
                   ufo::GeoVaLs &) const;

 private:
  void print(std::ostream &) const;
  F90getvalues keyGetValues_;
  ufo::Locations locs_;
  std::shared_ptr<const Geometry> geom_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
