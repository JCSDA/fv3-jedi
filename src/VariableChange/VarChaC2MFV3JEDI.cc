/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/VariableChange/VarChaC2MFV3JEDI.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "GeometryFV3JEDI.h"
#include "IncrementFV3JEDI.h"
#include "StateFV3JEDI.h"
#include "oops/util/Logger.h"
#include "VarChaC2MFV3JEDIFortran.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
VarChaC2MFV3JEDI::VarChaC2MFV3JEDI(const StateFV3JEDI & bg,
                                   const StateFV3JEDI & fg,
                                   const GeometryFV3JEDI & resol,
                                   const eckit::Configuration & conf):
    geom_(new GeometryFV3JEDI(resol))
{
    const eckit::Configuration * configc = &conf;
    fv3jedi_varcha_c2m_setup_f90(keyFtnConfig_, geom_->toFortran(),
                                 bg.toFortran(),
                                 fg.toFortran(),
                                 &configc);
    oops::Log::trace() << "VarChaC2MFV3JEDI created" << std::endl;
}
// -----------------------------------------------------------------------------
VarChaC2MFV3JEDI::~VarChaC2MFV3JEDI() {
    fv3jedi_varcha_c2m_delete_f90(keyFtnConfig_);
    oops::Log::trace() << "ChangeFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void VarChaC2MFV3JEDI::multiply(const IncrementFV3JEDI & dxa,
                                IncrementFV3JEDI & dxm) const {
  fv3jedi_varcha_c2m_multiply_f90(keyFtnConfig_, geom_->toFortran(),
                                  dxa.toFortran(),
                                  dxm.toFortran());
}
// -----------------------------------------------------------------------------
void VarChaC2MFV3JEDI::multiplyInverse(const IncrementFV3JEDI & dxm,
                                       IncrementFV3JEDI & dxa) const {
  fv3jedi_varcha_c2m_multiplyinverse_f90(keyFtnConfig_, geom_->toFortran(),
                                         dxm.toFortran(),
                                         dxa.toFortran());
}
// -----------------------------------------------------------------------------
void VarChaC2MFV3JEDI::multiplyAD(const IncrementFV3JEDI & dxm,
                                       IncrementFV3JEDI & dxa) const {
  fv3jedi_varcha_c2m_multiplyadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                         dxm.toFortran(),
                                         dxa.toFortran());
}
// -----------------------------------------------------------------------------
void VarChaC2MFV3JEDI::multiplyInverseAD(const IncrementFV3JEDI & dxa,
                                              IncrementFV3JEDI & dxm) const {
  fv3jedi_varcha_c2m_multiplyinverseadjoint_f90(keyFtnConfig_,
                                                geom_->toFortran(),
                                                dxm.toFortran(),
                                                dxa.toFortran());
}
// -----------------------------------------------------------------------------
void VarChaC2MFV3JEDI::print(std::ostream & os) const {
  os << "FV3JEDI change variable";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi

