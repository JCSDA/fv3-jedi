/*
 * (C) Copyright 2017-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/LinearVariableChange/Analysis2Model/LinVarChaA2M.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Utilities/Traits.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static LinearVariableChangeMaker<LinVarChaA2M> makerVarChaA2M_("Analysis2Model");
// -------------------------------------------------------------------------------------------------
LinVarChaA2M::LinVarChaA2M(const State & bg, const State & fg, const Geometry & resol,
                           const eckit::LocalConfiguration & conf)
  : LinearVariableChangeBase(), geom_(new Geometry(resol))
{
  util::Timer timer(classname(), "LinVarChaA2M");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  fv3jedi_linvarcha_a2m_create_f90(keyFtnConfig_, geom_->toFortran(), bg.toFortran(),
                                   fg.toFortran(), conf);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
LinVarChaA2M::~LinVarChaA2M() {
  util::Timer timer(classname(), "~LinVarChaA2M");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  fv3jedi_linvarcha_a2m_delete_f90(keyFtnConfig_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaA2M::multiply(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiply");
  oops::Log::trace() << classname() << " multiply starting" << std::endl;
  fv3jedi_linvarcha_a2m_multiply_f90(keyFtnConfig_, geom_->toFortran(), dxin.toFortran(),
                                     dxout.toFortran());
  oops::Log::trace() << classname() << " multiply done" << std::endl;
  dxout.validTime() = dxin.validTime();
}
// -------------------------------------------------------------------------------------------------
void LinVarChaA2M::multiplyInverse(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiplyInverse");
  oops::Log::trace() << classname() << " multiplyInverse starting" << std::endl;
  fv3jedi_linvarcha_a2m_multiplyinverse_f90(keyFtnConfig_, geom_->toFortran(), dxin.toFortran(),
                                            dxout.toFortran());
  oops::Log::trace() << classname() << " multiplyInverse done" << std::endl;
  dxout.validTime() = dxin.validTime();
}
// -------------------------------------------------------------------------------------------------
void LinVarChaA2M::multiplyAD(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiplyAD");
  oops::Log::trace() << classname() << " multiplyAD starting" << std::endl;
  fv3jedi_linvarcha_a2m_multiplyadjoint_f90(keyFtnConfig_, geom_->toFortran(), dxin.toFortran(),
                                            dxout.toFortran());
  oops::Log::trace() << classname() << " multiplyAD done" << std::endl;
  dxout.validTime() = dxin.validTime();
}
// -------------------------------------------------------------------------------------------------
void LinVarChaA2M::multiplyInverseAD(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiplyInverseAD");
  oops::Log::trace() << classname() << " multiplyInverseAD starting" << std::endl;
  fv3jedi_linvarcha_a2m_multiplyinverseadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                                   dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " multiplyInverseAD done" << std::endl;
  dxout.validTime() = dxin.validTime();
}
// -------------------------------------------------------------------------------------------------
void LinVarChaA2M::print(std::ostream & os) const {
  os << classname() << " variable change";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
