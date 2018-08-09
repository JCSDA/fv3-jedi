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
#include "oops/util/Logger.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
// static oops::VariableChangeMaker<FV3JEDITraits, ChangeVar>
// makerChVarFV3JEDI("FV3JEDICV");
// -----------------------------------------------------------------------------
ChangeVarFV3JEDI::ChangeVarFV3JEDI(const eckit::Configuration & conf) {
    const eckit::Configuration * configc = &conf;
    fv3jedi_changevar_setup_f90(keyFtnConfig_, &configc);
    oops::Log::trace() << "ChangeVarFV3JEDI created" << std::endl;
}
// -----------------------------------------------------------------------------
ChangeVarFV3JEDI::~ChangeVarFV3JEDI() {
    fv3jedi_changevar_delete_f90(keyFtnConfig_);
    oops::Log::trace() << "ChangeFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ChangeVarFV3JEDI::linearize(const StateFV3JEDI & other,
                                 const GeometryFV3JEDI & resol) {
    fv3jedi_changevar_linearize_f90(keyFtnConfig_,resol.toFortran(),
                                     other.fields().toFortran());
}
// -----------------------------------------------------------------------------
void ChangeVarFV3JEDI::transform(const IncrementFV3JEDI & dxa,
                                IncrementFV3JEDI & dxm) const {
  fv3jedi_changevar_transform_f90(keyFtnConfig_,dxa.fields().toFortran(),
                                  dxm.fields().toFortran());
}
// -----------------------------------------------------------------------------
void ChangeVarFV3JEDI::transformInverse(const IncrementFV3JEDI & dxm,
                                       IncrementFV3JEDI & dxa) const {
  dxa = dxm;
}
// -----------------------------------------------------------------------------
void ChangeVarFV3JEDI::transformAdjoint(const IncrementFV3JEDI & dxm,
                                       IncrementFV3JEDI & dxa) const {
  fv3jedi_changevar_transformadjoint_f90(keyFtnConfig_,dxm.fields().toFortran(),
                                         dxa.fields().toFortran());
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

