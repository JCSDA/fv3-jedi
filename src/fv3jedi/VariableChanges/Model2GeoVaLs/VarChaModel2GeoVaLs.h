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
#include "VarChaModel2GeoVaLs.interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class Geometry;
  class State;

// -------------------------------------------------------------------------------------------------

class VarChaModel2GeoVaLs: public util::Printable,
                           private util::ObjectCounter<VarChaModel2GeoVaLs> {
 public:
  static const std::string classname() {return "fv3jedi::VarChaModel2GeoVaLs";}

  explicit VarChaModel2GeoVaLs(const Geometry &, const eckit::Configuration &);
  ~VarChaModel2GeoVaLs();

  void changeVar(const State &, State &) const;
  void changeVarInverse(const State &, State &) const;

 private:
  F90m2g keyFtnConfig_;
  std::shared_ptr<const Geometry> geom_;
  void print(std::ostream &) const override;
};
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
