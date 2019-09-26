/*
 * (C) Copyright 2017-2019  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Analysis2Model/VarChaA2M.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/State/State.h"
#include "oops/util/Logger.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
VarChaA2M::VarChaA2M(const Geometry & resol, const eckit::Configuration & conf):
    geom_(new Geometry(resol))
{
  oops::Log::trace() << "VarChaA2M::VarChaA2M start" << std::endl;
  const eckit::Configuration * configc = &conf;
    fv3jedi_varcha_a2m_create_f90(keyFtnConfig_, geom_->toFortran(), &configc);
  oops::Log::trace() << "VarChaA2M::VarChaA2M done" << std::endl;
}
// -----------------------------------------------------------------------------
VarChaA2M::~VarChaA2M() {
  fv3jedi_varcha_a2m_delete_f90(keyFtnConfig_);
  oops::Log::trace() << "ChangeFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void VarChaA2M::changeVar(const State & xa,
                                       State & xm) const {
  oops::Log::trace() << "VarChaA2M::changeVar starting" << xm <<
                        std::endl;
  util::DateTime * vtime = &xm.validTime();
  fv3jedi_varcha_a2m_changevar_f90(keyFtnConfig_, geom_->toFortran(),
                                   xa.toFortran(), xm.toFortran(),
                                   &vtime);
  xm.validTime() = xa.validTime();
  oops::Log::trace() << "VarChaA2M::changeVar done" << xm << std::endl;
}
// -----------------------------------------------------------------------------
void VarChaA2M::changeVarInverse(const State & xm,
                                              State & xa) const {
  oops::Log::trace() << "VarChaA2M::changeVarInverse starting" <<xm <<
                        std::endl;
  util::DateTime * vtime = &xa.validTime();
  fv3jedi_varcha_a2m_changevarinverse_f90(keyFtnConfig_, geom_->toFortran(),
                                          xm.toFortran(), xa.toFortran(),
                                          &vtime);
  xa.validTime() = xm.validTime();
  oops::Log::trace() << "VarChaA2M::changeVarInverse done" << xm <<
                        std::endl;
}
// -----------------------------------------------------------------------------
void VarChaA2M::print(std::ostream & os) const {
  os << "VarChaA2M";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
