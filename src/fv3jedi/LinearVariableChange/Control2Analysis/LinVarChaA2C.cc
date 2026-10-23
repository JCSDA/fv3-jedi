/*
 * (C) Copyright 2017-2026 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/LinearVariableChange/Control2Analysis/LinVarChaA2C.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static LinearVariableChangeMaker<LinVarChaA2C> makerLinVarChaA2C_("Analysis2Control");
// -------------------------------------------------------------------------------------------------
LinVarChaA2C::LinVarChaA2C(const State & bg, const State & fg, const Geometry & resol,
                           const eckit::LocalConfiguration & conf)
  : LinearVariableChangeBase(), geom_(resol)
{
  util::Timer timer(classname(), "LinVarChaA2C");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  fv3jedi_linvarcha_c2a_create_f90(keyFtnConfig_, geom_.toFortran(), bg.toFortran(),
                                   fg.toFortran(), conf);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
LinVarChaA2C::~LinVarChaA2C() {
  util::Timer timer(classname(), "~LinVarChaA2C");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  fv3jedi_linvarcha_c2a_delete_f90(keyFtnConfig_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaA2C::multiply(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiply");
  oops::Log::trace() << classname() << " multiply starting" << std::endl;
  fv3jedi_linvarcha_c2a_multiplyinverse_f90(keyFtnConfig_, geom_.toFortran(),
                                            dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " multiply starting" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaA2C::multiplyAD(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiplyAD");
  oops::Log::trace() << classname() << " multiplyAD starting" << std::endl;
  fv3jedi_linvarcha_c2a_multiplyinverseadjoint_f90(keyFtnConfig_, geom_.toFortran(),
                                                   dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " multiplyAD starting" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaA2C::print(std::ostream & os) const {
  os << classname() << " variable change";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
