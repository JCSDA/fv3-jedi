/*
 * (C) Copyright 2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/Printable.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "LinVarChaModel2GeoVaLs.interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class Geometry;
  class State;
  class Increment;

// -------------------------------------------------------------------------------------------------

class LinVarChaModel2GeoVaLs: public util::Printable,
                              private util::ObjectCounter<LinVarChaModel2GeoVaLs> {
 public:
  static const std::string classname() {return "fv3jedi::LinVarChaModel2GeoVaLs";}

  explicit LinVarChaModel2GeoVaLs(const State &, const State &, const Geometry &,
                                  const eckit::Configuration &);
  ~LinVarChaModel2GeoVaLs();

  void multiply(const Increment &, Increment &) const;
  void multiplyInverse(const Increment &, Increment &) const;
  void multiplyAD(const Increment &, Increment &) const;
  void multiplyInverseAD(const Increment &, Increment &) const;

 private:
  std::shared_ptr<const Geometry> geom_;
  F90lm2g keyFtnConfig_;
  void print(std::ostream &) const override;
};
// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
