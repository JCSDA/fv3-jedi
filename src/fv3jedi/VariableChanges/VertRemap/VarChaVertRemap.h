/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "oops/util/Printable.h"
#include "VarChaVertRemap.interface.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class Geometry;
  class State;

// -------------------------------------------------------------------------------------------------
class VarChaVertRemap: public util::Printable {
 public:
  static const std::string classname() {return "fv3jedi::VarChaVertRemap";}
  explicit VarChaVertRemap(const Geometry &, const eckit::Configuration &);
  ~VarChaVertRemap();
  void changeVar(const State &, State &) const;
  void changeVarInverse(const State &, State &) const;

 private:
  std::shared_ptr<const Geometry> geom_;
  F90vc_VR keyFtn_;
  void print(std::ostream &) const override;
};
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
