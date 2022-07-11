/*
 * (C) Copyright 2017-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Utilities/Traits.h"
#include "fv3jedi/VariableChange/GEOSRstToBkg/VarChaGeosRst2Bkg.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static VariableChangeMaker<VarChaGeosRst2Bkg> makerVarChaA2M_("GeosRst2Bkg");
// -------------------------------------------------------------------------------------------------
VarChaGeosRst2Bkg::VarChaGeosRst2Bkg(const Geometry & resol, const eckit::LocalConfiguration & conf)
  : VariableChangeBase(), geom_(new Geometry(resol)) {
  util::Timer timer(classname(), "VarChaGeosRst2Bkg");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  fv3jedi_vc_geosrst2bkg_create_f90(keyFtnConfig_, geom_->toFortran(), conf);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
VarChaGeosRst2Bkg::~VarChaGeosRst2Bkg() {
  util::Timer timer(classname(), "~VarChaGeosRst2Bkg");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  fv3jedi_vc_geosrst2bkg_delete_f90(keyFtnConfig_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaGeosRst2Bkg::changeVar(const State & xin, State & xout) const {
  util::Timer timer(classname(), "changeVar");
  oops::Log::trace() << classname() << " changeVar starting" << std::endl;
  fv3jedi_vc_geosrst2bkg_changevar_f90(keyFtnConfig_, geom_->toFortran(), xin.toFortran(),
                                       xout.toFortran());
  xout.validTime() = xin.validTime();
  oops::Log::trace() << classname() << " changeVar done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaGeosRst2Bkg::changeVarInverse(const State & xin, State & xout) const {
  util::Timer timer(classname(), "changeVarInverse");
  oops::Log::trace() << classname() << " changeVarInverse starting" << std::endl;
  fv3jedi_vc_geosrst2bkg_changevarinverse_f90(keyFtnConfig_, geom_->toFortran(),
                                          xin.toFortran(), xout.toFortran());
  xout.validTime() = xin.validTime();
  oops::Log::trace() << classname() << " changeVarInverse done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaGeosRst2Bkg::print(std::ostream & os) const {
  os << classname() << " variable change";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
