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
#include "oops/util/Printable.h"
#include "VarChaGeosRst2Bkg.interface.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

class VarChaGeosRst2Bkg: public util::Printable {
 public:
  static const std::string classname() {return "fv3jedi::VarChaGeosRst2Bkg";}
  explicit VarChaGeosRst2Bkg(const Geometry &, const eckit::Configuration &);
  ~VarChaGeosRst2Bkg();
  void changeVar(const State &, State &) const;
  void changeVarInverse(const State &, State &) const;

 private:
  std::shared_ptr<const Geometry> geom_;
  F90vc_R2B keyFtnConfig_;
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
