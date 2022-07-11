/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_NORMGRADIENT_NORMGRADIENT_H_
#define FV3JEDI_NORMGRADIENT_NORMGRADIENT_H_

#include <memory>
#include <ostream>
#include <string>

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class Geometry;
  class Increment;
  class State;

// -----------------------------------------------------------------------------

class NormGradient : public util::Printable,
                    private util::ObjectCounter<NormGradient> {
 public:
  static const std::string classname() {return "fv3jedi::NormGradient";}

  NormGradient(const Geometry &, const State &, const eckit::Configuration &) {}
  virtual ~NormGradient() {}

  void apply(Increment &) const {}

// Private
 private:
  void print(std::ostream & os) const {os << " NormGradient: print not implemented yet.";}
};

// -----------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // FV3JEDI_NORMGRADIENT_NORMGRADIENT_H_
