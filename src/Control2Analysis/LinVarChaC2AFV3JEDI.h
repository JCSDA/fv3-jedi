/*
 * (C) Copyright 2017-2018  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SRC_CONTROL2ANALYSIS_LINVARCHAC2AFV3JEDI_H_
#define SRC_CONTROL2ANALYSIS_LINVARCHAC2AFV3JEDI_H_

#include <ostream>
#include <string>

#include "LinVarChaC2AFV3JEDI.interface.h"
#include "GeometryFV3JEDI.h"
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

class LinVarChaC2AFV3JEDI: public util::Printable {
 public:
  static const std::string classname() {return "fv3jedi::LinVarChaC2AFV3JEDI";}

  explicit LinVarChaC2AFV3JEDI(const StateFV3JEDI &, const StateFV3JEDI &,
                        const GeometryFV3JEDI &, const eckit::Configuration &);
  ~LinVarChaC2AFV3JEDI();

/// Perform linear multiplications
  void multiply(const IncrementFV3JEDI &, IncrementFV3JEDI &) const;
  void multiplyInverse(const IncrementFV3JEDI &, IncrementFV3JEDI &) const;
  void multiplyAD(const IncrementFV3JEDI &, IncrementFV3JEDI &) const;
  void multiplyInverseAD(const IncrementFV3JEDI &, IncrementFV3JEDI &) const;

 private:
  boost::shared_ptr<const GeometryFV3JEDI> geom_;
  F90lvcc2a keyFtnConfig_;
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // SRC_CONTROL2ANALYSIS_LINVARCHAC2AFV3JEDI_H_
