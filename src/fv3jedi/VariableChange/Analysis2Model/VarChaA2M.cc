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
#include "fv3jedi/VariableChange/Analysis2Model/VarChaA2M.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static VariableChangeMaker<VarChaA2M> makerVarChaA2M_("Analysis2Model");
// -------------------------------------------------------------------------------------------------
VarChaA2M::VarChaA2M(const Geometry & geom, const eckit::LocalConfiguration & conf)
  : VariableChangeBase(), geom_(new Geometry(geom)) {
  util::Timer timer(classname(), "VarChaA2M");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  fv3jedi_varcha_a2m_create_f90(keyFtnConfig_, geom_->toFortran(), conf);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
VarChaA2M::~VarChaA2M() {
  util::Timer timer(classname(), "~VarChaA2M");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  fv3jedi_varcha_a2m_delete_f90(keyFtnConfig_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaA2M::changeVar(const State & xin, State & xout) const {
  util::Timer timer(classname(), "changeVar");
  oops::Log::trace() << classname() << " changeVar done" << std::endl;
  fv3jedi_varcha_a2m_changevar_f90(keyFtnConfig_, geom_->toFortran(), xin.toFortran(),
                                   xout.toFortran());
  xout.validTime() = xin.validTime();
  oops::Log::trace() << classname() << " changeVar done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaA2M::changeVarInverse(const State & xin, State & xout) const {
  util::Timer timer(classname(), "changeVarInverse");
  oops::Log::trace() << classname() << " changeVarInverse starting" << std::endl;
  fv3jedi_varcha_a2m_changevarinverse_f90(keyFtnConfig_, geom_->toFortran(), xin.toFortran(),
                                          xout.toFortran());
  xout.validTime() = xin.validTime();
  oops::Log::trace() << classname() << " changeVarInverse done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void VarChaA2M::print(std::ostream & os) const {
  os << classname() << " variable change";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
