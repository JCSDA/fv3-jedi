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
#include "oops/base/VariableChangeBase.h"
#include "VarChaA2M.interface.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

class VarChaA2M: public oops::VariableChangeBase<Traits> {
 public:
  static const std::string classname() {return "fv3jedi::VarChaA2M";}
  VarChaA2M(const Geometry &, const eckit::Configuration &);
  ~VarChaA2M();
  void changeVar(const State &, State &) const override;
  void changeVarInverse(const State &, State &) const override;

 private:
  std::shared_ptr<const Geometry> geom_;
  F90vc_A2M keyFtnConfig_;
  void print(std::ostream &) const override;
  eckit::LocalConfiguration conf_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
