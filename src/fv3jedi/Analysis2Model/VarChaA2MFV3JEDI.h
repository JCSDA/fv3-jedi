/*
 * (C) Copyright 2017-2019  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_ANALYSIS2MODEL_VARCHAA2MFV3JEDI_H_
#define FV3JEDI_ANALYSIS2MODEL_VARCHAA2MFV3JEDI_H_

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/GeometryFV3JEDI.h"
#include "oops/util/Printable.h"
#include "VarChaA2MFV3JEDI.interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class GeometryFV3JEDI;
  class StateFV3JEDI;

// -----------------------------------------------------------------------------
/// FV3JEDI nonlinear change of variable

class VarChaA2MFV3JEDI: public util::Printable {
 public:
  static const std::string classname() {return "fv3jedi::VarChaA2MFV3JEDI";}

  explicit VarChaA2MFV3JEDI(const GeometryFV3JEDI &,
                            const eckit::Configuration &);
  ~VarChaA2MFV3JEDI();

  void changeVar(const StateFV3JEDI &, StateFV3JEDI &) const;
  void changeVarInverse(const StateFV3JEDI &, StateFV3JEDI &) const;

 private:
  boost::shared_ptr<const GeometryFV3JEDI> geom_;
  F90vca2m keyFtnConfig_;
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------
}  // namespace fv3jedi

#endif  // FV3JEDI_ANALYSIS2MODEL_VARCHAA2MFV3JEDI_H_
