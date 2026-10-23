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
#include "fv3jedi/VariableChange/Analysis2Control/VarChaA2C.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static VariableChangeMaker<VarChaA2C> makerVarChaA2C_("Analysis2Control");
// -------------------------------------------------------------------------------------------------
VarChaA2C::VarChaA2C(const Geometry & resol, const eckit::LocalConfiguration & conf)
  : VariableChangeBase(), geom_(resol) {
  util::Timer timer(classname(), "VarChaA2C");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  fv3jedi_varcha_c2a_create_f90(keyFtnConfig_, geom_.toFortran(), conf);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
VarChaA2C::~VarChaA2C() {
  util::Timer timer(classname(), "~VarChaA2C");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  fv3jedi_varcha_c2a_delete_f90(keyFtnConfig_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaA2C::changeVar(const State & xin, State & xout) const {
  util::Timer timer(classname(), "changeVar");
  oops::Log::trace() << classname() << " changeVar starting" << std::endl;
  fv3jedi_varcha_c2a_changevarinverse_f90(keyFtnConfig_, geom_.toFortran(), xin.toFortran(),
                                          xout.toFortran());
  xout.validTime() = xin.validTime();
  oops::Log::trace() << classname() << " changeVar done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaA2C::print(std::ostream & os) const {
  os << classname() << " variable change";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
