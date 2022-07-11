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
#include "fv3jedi/Utilities/Traits.h"
#include "fv3jedi/VariableChange/Base/VariableChangeBase.h"
#include "VarChaColdStartWinds.interface.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

class VarChaColdStartWinds: public VariableChangeBase {
 public:
  static const std::string classname() {return "fv3jedi::VarChaColdStartWinds";}
  VarChaColdStartWinds(const Geometry &, const eckit::LocalConfiguration &);
  ~VarChaColdStartWinds();
  void changeVar(const State &, State &) const override;
  void changeVarInverse(const State &, State &) const override;

 private:
  F90vc_CSW keyFtn_;
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
