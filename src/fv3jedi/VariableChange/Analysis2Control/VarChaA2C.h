/*
 * (C) Copyright 2017-2020  UCAR.
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
#include "fv3jedi/VariableChange/Base/VariableChangeBase.h"
#include "fv3jedi/VariableChange/Control2Analysis/VarChaC2A.interface.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

class VarChaA2C: public VariableChangeBase {
 public:
  static const std::string classname() {return "fv3jedi::VarChaA2C";}
  VarChaA2C(const Geometry &, const eckit::LocalConfiguration &);
  ~VarChaA2C();
  void changeVar(const State &, State &) const override;

 private:
  const Geometry & geom_;
  F90vc_C2A keyFtnConfig_;
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
