/*
 * (C) Copyright 2017-2019  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_VARIABLECHANGES_GEOSRSTTOBKG_VARCHAGEOSRST2BKG_H_
#define FV3JEDI_VARIABLECHANGES_GEOSRSTTOBKG_VARCHAGEOSRST2BKG_H_

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "oops/util/Printable.h"
#include "VarChaGeosRst2Bkg.interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class Geometry;
  class State;

// -----------------------------------------------------------------------------
/// FV3JEDI nonlinear change of variable

class VarChaGeosRst2Bkg: public util::Printable {
 public:
  static const std::string classname() {return "fv3jedi::VarChaGeosRst2Bkg";}

  explicit VarChaGeosRst2Bkg(const Geometry &,
                            const eckit::Configuration &);
  ~VarChaGeosRst2Bkg();

  void changeVar(const State &, State &) const;
  void changeVarInverse(const State &, State &) const;

 private:
  std::shared_ptr<const Geometry> geom_;
  F90vcd2a keyFtnConfig_;
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------
}  // namespace fv3jedi

#endif  // FV3JEDI_VARIABLECHANGES_GEOSRSTTOBKG_VARCHAGEOSRST2BKG_H_
