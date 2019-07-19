/*
 * (C) Copyright 2017-2019  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_ANALYSIS2MODEL_VARCHAA2M_H_
#define FV3JEDI_ANALYSIS2MODEL_VARCHAA2M_H_

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "oops/util/Printable.h"
#include "VarChaA2M.interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class Geometry;
  class State;

// -----------------------------------------------------------------------------
/// FV3JEDI nonlinear change of variable

class VarChaA2M: public util::Printable {
 public:
  static const std::string classname() {return "fv3jedi::VarChaA2M";}

  explicit VarChaA2M(const Geometry &,
                            const eckit::Configuration &);
  ~VarChaA2M();

  void changeVar(const State &, State &) const;
  void changeVarInverse(const State &, State &) const;

 private:
  boost::shared_ptr<const Geometry> geom_;
  F90vca2m keyFtnConfig_;
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------
}  // namespace fv3jedi

#endif  // FV3JEDI_ANALYSIS2MODEL_VARCHAA2M_H_
