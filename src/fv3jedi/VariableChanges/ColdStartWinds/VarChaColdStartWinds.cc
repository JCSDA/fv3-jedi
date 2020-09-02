/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/VariableChanges/ColdStartWinds/VarChaColdStartWinds.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/State/State.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
VarChaColdStartWinds::VarChaColdStartWinds(const Geometry & resol,
                                           const eckit::Configuration & conf) {
  util::Timer timer(classname(), "VarChaColdStartWinds");
  oops::Log::trace() << "VarChaColdStartWinds::VarChaColdStartWinds start" << std::endl;
  const eckit::Configuration * configc = &conf;
  fv3jedi_vc_coldstartwinds_create_f90(keyFtn_, resol.toFortran(), &configc);
  oops::Log::trace() << "VarChaColdStartWinds::VarChaColdStartWinds done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
VarChaColdStartWinds::~VarChaColdStartWinds() {
  util::Timer timer(classname(), "~VarChaColdStartWinds");
  oops::Log::trace() << "VarChaColdStartWinds::~VarChaColdStartWinds start" << std::endl;
  fv3jedi_vc_coldstartwinds_delete_f90(keyFtn_);
  oops::Log::trace() << "VarChaColdStartWinds::~VarChaColdStartWinds done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaColdStartWinds::changeVar(const State & xin, State & xout) const {
  util::Timer timer(classname(), "changeVar");
  oops::Log::trace() << "VarChaColdStartWinds::changeVar starting" << xin << std::endl;
  fv3jedi_vc_coldstartwinds_changevar_f90(keyFtn_, xin.toFortran(), xout.toFortran());
  xout.validTime() = xin.validTime();
  oops::Log::trace() << "VarChaColdStartWinds::changeVar done" << xout << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaColdStartWinds::changeVarInverse(const State & xin, State & xout) const {
  util::Timer timer(classname(), "changeVarInverse");
  oops::Log::trace() << "VarChaColdStartWinds::changeVarInverse starting" << xin << std::endl;
  xout = xin;  // No inverse required
  xout.validTime() = xin.validTime();
  oops::Log::trace() << "VarChaColdStartWinds::changeVarInverse done" << xout << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaColdStartWinds::print(std::ostream & os) const {
  os << "VarChaColdStartWinds";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
