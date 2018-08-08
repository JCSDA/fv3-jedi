/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/ChangeVarFV3JEDI.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "src/GeometryFV3JEDI.h"
#include "src/IncrementFV3JEDI.h"
#include "src/StateFV3JEDI.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
// static oops::VariableChangeMaker<FV3JEDITraits, ChangeVar>
// makerChVarFV3JEDI("FV3JEDICV");
// -----------------------------------------------------------------------------
ChangeVarFV3JEDI::ChangeVarFV3JEDI(const eckit::Configuration &) {}
// -----------------------------------------------------------------------------
ChangeVarFV3JEDI::~ChangeVarFV3JEDI() {}
// -----------------------------------------------------------------------------
void ChangeVarFV3JEDI::linearize(const StateFV3JEDI &,
                                 const GeometryFV3JEDI &) {}
// -----------------------------------------------------------------------------
void ChangeVarFV3JEDI::transform(const IncrementFV3JEDI & dxa,
                                IncrementFV3JEDI & dxm) const {
  dxm = dxa;
}
// -----------------------------------------------------------------------------
void ChangeVarFV3JEDI::transformInverse(const IncrementFV3JEDI & dxm,
                                       IncrementFV3JEDI & dxa) const {
  dxa = dxm;
}
// -----------------------------------------------------------------------------
void ChangeVarFV3JEDI::transformAdjoint(const IncrementFV3JEDI & dxm,
                                       IncrementFV3JEDI & dxa) const {
  dxa = dxm;
}
// -----------------------------------------------------------------------------
void ChangeVarFV3JEDI::transformInverseAdjoint(const IncrementFV3JEDI & dxa,
                                              IncrementFV3JEDI & dxm) const {
  dxm = dxa;
}
// -----------------------------------------------------------------------------
void ChangeVarFV3JEDI::print(std::ostream & os) const {
  os << "FV3JEDI change variable";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi

