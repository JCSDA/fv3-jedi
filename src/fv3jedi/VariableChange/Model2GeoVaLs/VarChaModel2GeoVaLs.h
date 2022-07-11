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
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Utilities/Traits.h"
#include "fv3jedi/VariableChange/Base/VariableChangeBase.h"
#include "VarChaModel2GeoVaLs.interface.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

class VarChaModel2GeoVaLs: public VariableChangeBase,
                           private util::ObjectCounter<VarChaModel2GeoVaLs> {
 public:
  static const std::string classname() {return "fv3jedi::VarChaModel2GeoVaLs";}
  VarChaModel2GeoVaLs(const Geometry &, const eckit::LocalConfiguration &);
  ~VarChaModel2GeoVaLs();
  void changeVar(const State &, State &) const override;
  void changeVarInverse(const State &, State &) const override;

 private:
  F90vc_M2G keyFtnConfig_;
  std::shared_ptr<const Geometry> geom_;
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
