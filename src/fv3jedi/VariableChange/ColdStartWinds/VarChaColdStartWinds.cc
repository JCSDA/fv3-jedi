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
#include "fv3jedi/VariableChange/ColdStartWinds/VarChaColdStartWinds.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static VariableChangeMaker<VarChaColdStartWinds> makerVarChaA2M_("ColdStartWinds");
// -------------------------------------------------------------------------------------------------
VarChaColdStartWinds::VarChaColdStartWinds(const Geometry & resol,
                                           const eckit::LocalConfiguration & conf)
  : VariableChangeBase() {
  util::Timer timer(classname(), "VarChaColdStartWinds");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  fv3jedi_vc_coldstartwinds_create_f90(keyFtn_, resol.toFortran(), conf);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
VarChaColdStartWinds::~VarChaColdStartWinds() {
  util::Timer timer(classname(), "~VarChaColdStartWinds");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  fv3jedi_vc_coldstartwinds_delete_f90(keyFtn_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaColdStartWinds::changeVar(const State & xin, State & xout) const {
  util::Timer timer(classname(), "changeVar");
  oops::Log::trace() << classname() << " changeVar starting" << std::endl;
  fv3jedi_vc_coldstartwinds_changevar_f90(keyFtn_, xin.toFortran(), xout.toFortran());
  xout.validTime() = xin.validTime();
  oops::Log::trace() << classname() << " changeVar done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaColdStartWinds::changeVarInverse(const State & xin, State & xout) const {
  util::Timer timer(classname(), "changeVarInverse");
  oops::Log::trace() << classname() << " changeVarInverse starting" << std::endl;
  xout = xin;  // No inverse required
  xout.validTime() = xin.validTime();
  oops::Log::trace() << classname() << " changeVarInverse done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaColdStartWinds::print(std::ostream & os) const {
  os << classname() << " variable change";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
