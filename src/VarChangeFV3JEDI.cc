/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/VarChangeFV3JEDI.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "src/GeometryFV3JEDI.h"
#include "src/IncrementFV3JEDI.h"
#include "src/StateFV3JEDI.h"
#include "oops/util/Logger.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
VarChangeFV3JEDI::VarChangeFV3JEDI(const eckit::Configuration & conf) {
    const eckit::Configuration * configc = &conf;
    fv3jedi_varchange_setup_f90(keyFtnConfig_, &configc);
    oops::Log::trace() << "VarChangeFV3JEDI created" << std::endl;
}
// -----------------------------------------------------------------------------
VarChangeFV3JEDI::~VarChangeFV3JEDI() {
    fv3jedi_varchange_delete_f90(keyFtnConfig_);
    oops::Log::trace() << "ChangeFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void VarChangeFV3JEDI::linearize(const StateFV3JEDI & other,
                                 const GeometryFV3JEDI & resol) {
    fv3jedi_varchange_linearize_f90(keyFtnConfig_,resol.toFortran(),
                                     other.fields().toFortran());
}
// -----------------------------------------------------------------------------
void VarChangeFV3JEDI::multiply(const IncrementFV3JEDI & dxa,
                                IncrementFV3JEDI & dxm) const {
  fv3jedi_varchange_multiply_f90(keyFtnConfig_,dxa.fields().toFortran(),
                                  dxm.fields().toFortran());
}
// -----------------------------------------------------------------------------
void VarChangeFV3JEDI::multiplyInverse(const IncrementFV3JEDI & dxm,
                                       IncrementFV3JEDI & dxa) const {
  dxa = dxm;
}
// -----------------------------------------------------------------------------
void VarChangeFV3JEDI::multiplyAD(const IncrementFV3JEDI & dxm,
                                       IncrementFV3JEDI & dxa) const {
  fv3jedi_varchange_multiplyadjoint_f90(keyFtnConfig_,dxm.fields().toFortran(),
                                         dxa.fields().toFortran());
}
// -----------------------------------------------------------------------------
void VarChangeFV3JEDI::multiplyInverseAD(const IncrementFV3JEDI & dxa,
                                              IncrementFV3JEDI & dxm) const {
  dxm = dxa;
}
// -----------------------------------------------------------------------------
void VarChangeFV3JEDI::print(std::ostream & os) const {
  os << "FV3JEDI change variable";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi

