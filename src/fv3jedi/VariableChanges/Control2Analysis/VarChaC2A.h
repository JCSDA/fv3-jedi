/*
 * (C) Copyright 2017-2019  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_VARIABLECHANGES_CONTROL2ANALYSIS_VARCHAC2A_H_
#define FV3JEDI_VARIABLECHANGES_CONTROL2ANALYSIS_VARCHAC2A_H_

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "oops/util/Printable.h"
#include "VarChaC2A.interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class Geometry;
  class State;

// -----------------------------------------------------------------------------
/// FV3JEDI nonlinear change of variable

class VarChaC2A: public util::Printable {
 public:
  static const std::string classname() {return "fv3jedi::VarChaC2A";}

  explicit VarChaC2A(const Geometry &,
                            const eckit::Configuration &);
  ~VarChaC2A();

  void changeVar(const State &, State &) const;
  void changeVarInverse(const State &, State &) const;

 private:
  boost::shared_ptr<const Geometry> geom_;
  F90vcc2a keyFtnConfig_;
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------
}  // namespace fv3jedi

#endif  // FV3JEDI_VARIABLECHANGES_CONTROL2ANALYSIS_VARCHAC2A_H_
