/*
 * (C) Copyright 2017-2020  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/interface/LinearVariableChange.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Utilities/Traits.h"
#include "fv3jedi/Utilities/Utilities.h"
#include "fv3jedi/VariableChanges/Control2Analysis/LinVarChaC2A.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static oops::LinearVariableChangeMaker<Traits, oops::LinearVariableChange<Traits, LinVarChaC2A> >
       makerLinVarChaC2A_("Control2Analysis");
// -------------------------------------------------------------------------------------------------
LinVarChaC2A::LinVarChaC2A(const State & bg, const State & fg, const Geometry & resol,
                           const eckit::Configuration & conf):
  geom_(new Geometry(resol))
{
  util::Timer timer(classname(), "LinVarChaC2A");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  const eckit::Configuration * configc = &conf;
  fv3jedi_linvarcha_c2a_create_f90(keyFtnConfig_, geom_->toFortran(), bg.toFortran(),
                                   fg.toFortran(), &configc);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
LinVarChaC2A::~LinVarChaC2A() {
  util::Timer timer(classname(), "~LinVarChaC2A");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  fv3jedi_linvarcha_c2a_delete_f90(keyFtnConfig_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaC2A::multiply(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiply");
  oops::Log::trace() << classname() << " multiply starting" << std::endl;
  fv3jedi_linvarcha_c2a_multiply_f90(keyFtnConfig_, geom_->toFortran(),
                                     dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " multiply starting" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaC2A::multiplyInverse(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiplyInverse");
  oops::Log::trace() << classname() << " multiplyInverse starting" << std::endl;
  fv3jedi_linvarcha_c2a_multiplyinverse_f90(keyFtnConfig_, geom_->toFortran(),
                                            dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " multiplyInverse starting" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaC2A::multiplyAD(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiplyAD");
  oops::Log::trace() << classname() << " multiplyAD starting" << std::endl;
  fv3jedi_linvarcha_c2a_multiplyadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                            dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " multiplyAD starting" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaC2A::multiplyInverseAD(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiplyInverseAD");
  oops::Log::trace() << classname() << " multiplyInverseAD starting" << std::endl;
  fv3jedi_linvarcha_c2a_multiplyinverseadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                                   dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " multiplyInverseAD starting" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaC2A::print(std::ostream & os) const {
  os << classname() << " variable change";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
