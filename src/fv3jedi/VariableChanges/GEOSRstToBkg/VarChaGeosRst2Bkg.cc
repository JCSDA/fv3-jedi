/*
 * (C) Copyright 2017-2019  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/VariableChanges/GEOSRstToBkg/VarChaGeosRst2Bkg.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/State/State.h"
#include "oops/util/Logger.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
VarChaGeosRst2Bkg::VarChaGeosRst2Bkg(const Geometry & resol, const eckit::Configuration & conf):
    geom_(new Geometry(resol))
{
  oops::Log::trace() << "VarChaGeosRst2Bkg::VarChaGeosRst2Bkg start" << std::endl;
  const eckit::Configuration * configc = &conf;
    fv3jedi_vc_geosrst2bkg_create_f90(keyFtnConfig_, geom_->toFortran(), &configc);
  oops::Log::trace() << "VarChaGeosRst2Bkg::VarChaGeosRst2Bkg done" << std::endl;
}
// -----------------------------------------------------------------------------
VarChaGeosRst2Bkg::~VarChaGeosRst2Bkg() {
  fv3jedi_vc_geosrst2bkg_delete_f90(keyFtnConfig_);
  oops::Log::trace() << "ChangeFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void VarChaGeosRst2Bkg::changeVar(const State & xd, State & xa) const {
  oops::Log::trace() << "VarChaGeosRst2Bkg::changeVar starting" << xa << std::endl;
  fv3jedi_vc_geosrst2bkg_changevar_f90(keyFtnConfig_, geom_->toFortran(),
                                   xd.toFortran(), xa.toFortran());
  xa.validTime() = xd.validTime();
  oops::Log::trace() << "VarChaGeosRst2Bkg::changeVar done" << xa << std::endl;
}
// -----------------------------------------------------------------------------
void VarChaGeosRst2Bkg::changeVarInverse(const State & xa, State & xd) const {
  oops::Log::trace() << "VarChaGeosRst2Bkg::changeVarInverse starting" <<xa << std::endl;
  fv3jedi_vc_geosrst2bkg_changevarinverse_f90(keyFtnConfig_, geom_->toFortran(),
                                          xa.toFortran(), xd.toFortran());
  xd.validTime() = xa.validTime();
  oops::Log::trace() << "VarChaGeosRst2Bkg::changeVarInverse done" << xa <<
                        std::endl;
}
// -----------------------------------------------------------------------------
void VarChaGeosRst2Bkg::print(std::ostream & os) const {
  os << "VarChaGeosRst2Bkg";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
