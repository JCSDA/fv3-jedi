/*
 * (C) Copyright 2017-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "fv3jedi/LinearVariableChange/Base/LinearVariableChangeBase.h"
#include "fv3jedi/LinearVariableChange/Model2GeoVaLs/LinVarChaModel2GeoVaLs.interface.h"
#include "fv3jedi/Utilities/Traits.h"

namespace fv3jedi {
  class Geometry;
  class Increment;
  class State;

// -------------------------------------------------------------------------------------------------

class LinVarChaModel2GeoVaLs : public LinearVariableChangeBase {
 public:
  static const std::string classname() {return "fv3jedi::LinVarChaModel2GeoVaLs";}
  explicit LinVarChaModel2GeoVaLs(const State &, const State &, const Geometry &,
                                  const eckit::LocalConfiguration &);
  ~LinVarChaModel2GeoVaLs();
  void multiply(const Increment &, Increment &) const override;
  void multiplyInverse(const Increment &, Increment &) const override;
  void multiplyAD(const Increment &, Increment &) const override;
  void multiplyInverseAD(const Increment &, Increment &) const override;

 private:
  std::shared_ptr<const Geometry> geom_;
  F90lvc_M2G keyFtnConfig_;
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
