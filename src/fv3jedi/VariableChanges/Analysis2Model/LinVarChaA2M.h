/*
 * (C) Copyright 2017-2018  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_VARIABLECHANGES_ANALYSIS2MODEL_LINVARCHAA2M_H_
#define FV3JEDI_VARIABLECHANGES_ANALYSIS2MODEL_LINVARCHAA2M_H_

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "oops/util/Printable.h"
#include "LinVarChaA2M.interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class Geometry;
  class State;
  class Increment;

// -----------------------------------------------------------------------------
/// FV3JEDI linear change of variable

class LinVarChaA2M: public util::Printable {
 public:
  static const std::string classname() {return "fv3jedi::LinVarChaA2M";}

  explicit LinVarChaA2M(const State &, const State &,
                        const Geometry &, const eckit::Configuration &);
  ~LinVarChaA2M();

/// Perform linear multiplications
  void multiply(const Increment &, Increment &) const;
  void multiplyInverse(const Increment &, Increment &) const;
  void multiplyAD(const Increment &, Increment &) const;
  void multiplyInverseAD(const Increment &, Increment &) const;

 private:
  boost::shared_ptr<const Geometry> geom_;
  F90lvca2m keyFtnConfig_;
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_VARIABLECHANGES_ANALYSIS2MODEL_LINVARCHAA2M_H_
