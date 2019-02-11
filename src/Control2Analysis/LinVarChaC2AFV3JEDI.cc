/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Control2Analysis/LinVarChaC2AFV3JEDI.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "GeometryFV3JEDI.h"
#include "IncrementFV3JEDI.h"
#include "StateFV3JEDI.h"
#include "oops/util/Logger.h"
#include "LinVarChaC2AFV3JEDI.interface.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
LinVarChaC2AFV3JEDI::LinVarChaC2AFV3JEDI(const StateFV3JEDI & bg,
                                   const StateFV3JEDI & fg,
                                   const GeometryFV3JEDI & resol,
                                   const eckit::Configuration & conf):
    geom_(new GeometryFV3JEDI(resol))
{
    const eckit::Configuration * configc = &conf;
    fv3jedi_linvarcha_c2a_setup_f90(keyFtnConfig_, geom_->toFortran(),
                                 bg.toFortran(),
                                 fg.toFortran(),
                                 &configc);
    oops::Log::trace() << "LinVarChaC2AFV3JEDI created" << std::endl;
}
// -----------------------------------------------------------------------------
LinVarChaC2AFV3JEDI::~LinVarChaC2AFV3JEDI() {
    fv3jedi_linvarcha_c2a_delete_f90(keyFtnConfig_);
    oops::Log::trace() << "ChangeFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void LinVarChaC2AFV3JEDI::multiply(const IncrementFV3JEDI & dxc,
                                         IncrementFV3JEDI & dxa) const {
  fv3jedi_linvarcha_c2a_multiply_f90(keyFtnConfig_, geom_->toFortran(),
                                     dxc.toFortran(), dxa.toFortran());
}
// -----------------------------------------------------------------------------
void LinVarChaC2AFV3JEDI::multiplyInverse(const IncrementFV3JEDI & dxa,
                                                IncrementFV3JEDI & dxc) const {
  fv3jedi_linvarcha_c2a_multiplyinverse_f90(keyFtnConfig_, geom_->toFortran(),
                                            dxa.toFortran(), dxc.toFortran());
}
// -----------------------------------------------------------------------------
void LinVarChaC2AFV3JEDI::multiplyAD(const IncrementFV3JEDI & dxa,
                                           IncrementFV3JEDI & dxc) const {
  fv3jedi_linvarcha_c2a_multiplyadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                            dxa.toFortran(), dxc.toFortran());
}
// -----------------------------------------------------------------------------
void LinVarChaC2AFV3JEDI::multiplyInverseAD(const IncrementFV3JEDI & dxc,
                                                  IncrementFV3JEDI & dxa) const {
  fv3jedi_linvarcha_c2a_multiplyinverseadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                                   dxc.toFortran(), dxa.toFortran());
}
// -----------------------------------------------------------------------------
void LinVarChaC2AFV3JEDI::print(std::ostream & os) const {
  os << "FV3JEDI change variable";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
