/*
 * (C) Copyright 2017-2019  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "src/Analysis2Model/VarChaA2MFV3JEDI.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "GeometryFV3JEDI.h"
#include "StateFV3JEDI.h"
#include "oops/util/Logger.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
VarChaA2MFV3JEDI::VarChaA2MFV3JEDI(const GeometryFV3JEDI & resol,
                                   const eckit::Configuration & conf):
    geom_(new GeometryFV3JEDI(resol))
{
  oops::Log::trace() << "VarChaA2MFV3JEDI::VarChaA2MFV3JEDI start" << std::endl;
  const eckit::Configuration * configc = &conf;
    fv3jedi_varcha_a2m_create_f90(keyFtnConfig_, geom_->toFortran(), &configc);
  oops::Log::trace() << "VarChaA2MFV3JEDI::VarChaA2MFV3JEDI done" << std::endl;
}
// -----------------------------------------------------------------------------
VarChaA2MFV3JEDI::~VarChaA2MFV3JEDI() {
  fv3jedi_varcha_a2m_delete_f90(keyFtnConfig_);
  oops::Log::trace() << "ChangeFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void VarChaA2MFV3JEDI::changeVar(const StateFV3JEDI & xa,
                                       StateFV3JEDI & xm) const {
  oops::Log::trace() << "VarChaA2MFV3JEDI::changeVar starting" << xm <<
                        std::endl;
  util::DateTime * vtime = &xm.validTime();
  fv3jedi_varcha_a2m_changevar_f90(keyFtnConfig_, geom_->toFortran(),
                                   xa.toFortran(), xm.toFortran(),
                                   &vtime);
  xm.validTime() = xa.validTime();
  oops::Log::trace() << "VarChaA2MFV3JEDI::changeVar done" << xm << std::endl;
}
// -----------------------------------------------------------------------------
void VarChaA2MFV3JEDI::changeVarInverse(const StateFV3JEDI & xm,
                                              StateFV3JEDI & xa) const {
  oops::Log::trace() << "VarChaA2MFV3JEDI::changeVarInverse starting" <<xm <<
                        std::endl;
  util::DateTime * vtime = &xa.validTime();
  fv3jedi_varcha_a2m_changevarinverse_f90(keyFtnConfig_, geom_->toFortran(),
                                          xm.toFortran(), xa.toFortran(),
                                          &vtime);
  xa.validTime() = xm.validTime();
  oops::Log::trace() << "VarChaA2MFV3JEDI::changeVarInverse done" << xm <<
                        std::endl;
}
// -----------------------------------------------------------------------------
void VarChaA2MFV3JEDI::print(std::ostream & os) const {
  os << "VarChaA2MFV3JEDI";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
