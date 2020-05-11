/*
 * (C) Copyright 2017-2019  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/VariableChanges/Control2Analysis/VarChaC2A.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/State/State.h"
#include "oops/util/Logger.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
VarChaC2A::VarChaC2A(const Geometry & resol, const eckit::Configuration & conf):
    geom_(new Geometry(resol))
{
  oops::Log::trace() << "VarChaC2A::VarChaC2A start" << std::endl;
  const eckit::Configuration * configc = &conf;
    fv3jedi_varcha_c2a_create_f90(keyFtnConfig_, geom_->toFortran(), &configc);
  oops::Log::trace() << "VarChaC2A::VarChaC2A done" << std::endl;
}
// -----------------------------------------------------------------------------
VarChaC2A::~VarChaC2A() {
  fv3jedi_varcha_c2a_delete_f90(keyFtnConfig_);
  oops::Log::trace() << "VarChaC2A destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void VarChaC2A::changeVar(const State & xa, State & xm) const {
  oops::Log::trace() << "VarChaC2A::changeVar starting" << xm <<
                        std::endl;
  fv3jedi_varcha_c2a_changevar_f90(keyFtnConfig_, geom_->toFortran(),
                                   xa.toFortran(), xm.toFortran());
  xm.validTime() = xa.validTime();
  oops::Log::trace() << "VarChaC2A::changeVar done" << xm << std::endl;
}
// -----------------------------------------------------------------------------
void VarChaC2A::changeVarInverse(const State & xm, State & xa) const {
  oops::Log::trace() << "VarChaC2A::changeVarInverse starting" <<xm <<
                        std::endl;
  fv3jedi_varcha_c2a_changevarinverse_f90(keyFtnConfig_, geom_->toFortran(),
                                          xm.toFortran(), xa.toFortran());
  xa.validTime() = xm.validTime();
  oops::Log::trace() << "VarChaC2A::changeVarInverse done" << xm <<
                        std::endl;
}
// -----------------------------------------------------------------------------
void VarChaC2A::print(std::ostream & os) const {
  os << "VarChaC2A";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
