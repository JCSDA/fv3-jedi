/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_MODEL_CHANGEVAR_H_
#define FV3JEDI_MODEL_CHANGEVAR_H_

#include <ostream>
#include <string>

#include "Fortran.h"
#include "eckit/config/Configuration.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class GeometryFV3JEDI;
  class StateFV3JEDI;
  class IncrementFV3JEDI;

// -----------------------------------------------------------------------------
/// FV3JEDI linear change of variable

class ChangeVarFV3JEDI: public util::Printable {
 public:
  static const std::string classname() {return "fv3jedi::ChangeVar";}

  explicit ChangeVarFV3JEDI(const eckit::Configuration &);
  ~ChangeVarFV3JEDI();

/// Set linearisation state
  void linearize(const StateFV3JEDI &, const GeometryFV3JEDI &);

/// Perform linear transforms
  void transform(const IncrementFV3JEDI &, IncrementFV3JEDI &) const;
  void transformInverse(const IncrementFV3JEDI &, IncrementFV3JEDI &) const;
  void transformAdjoint(const IncrementFV3JEDI &, IncrementFV3JEDI &) const;
  void transformInverseAdjoint(const IncrementFV3JEDI &,
                                     IncrementFV3JEDI &) const;

 private:
  F90cvar keyFtnConfig_;
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_MODEL_CHANGEVAR_H_
