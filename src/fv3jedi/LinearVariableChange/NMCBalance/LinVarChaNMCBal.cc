/*
 * (C) Copyright 2017-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/LinearVariableChange/NMCBalance/LinVarChaNMCBal.h"

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
static LinearVariableChangeMaker<LinVarChaNMCBal> makerLinVarChaNMCBal_("NMCBalance");
// -------------------------------------------------------------------------------------------------
LinVarChaNMCBal::LinVarChaNMCBal(const State & bg, const State & fg, const Geometry & resol,
                                 const eckit::LocalConfiguration & conf)
  : LinearVariableChangeBase(), geom_(new Geometry(resol))
{
  util::Timer timer(classname(), "LinVarChaNMCBal");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;
  fv3jedi_linvarcha_nmcbal_create_f90(keyFtnConfig_, geom_->toFortran(), bg.toFortran(),
                                      fg.toFortran(), conf);
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
LinVarChaNMCBal::~LinVarChaNMCBal() {
  util::Timer timer(classname(), "~LinVarChaNMCBal");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  fv3jedi_linvarcha_nmcbal_delete_f90(keyFtnConfig_);
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaNMCBal::multiply(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiply");
  oops::Log::trace() << classname() << " multiply starting" << std::endl;
  fv3jedi_linvarcha_nmcbal_multiply_f90(keyFtnConfig_, geom_->toFortran(),
                                        dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " multiply done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaNMCBal::multiplyInverse(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiplyInverse");
  oops::Log::trace() << classname() << " multiplyInverse starting" << std::endl;
  fv3jedi_linvarcha_nmcbal_multiplyinverse_f90(keyFtnConfig_, geom_->toFortran(),
                                               dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " multiplyInverse done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaNMCBal::multiplyAD(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiplyAD");
  oops::Log::trace() << classname() << " multiplyAD starting" << std::endl;
  fv3jedi_linvarcha_nmcbal_multiplyadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                               dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " multiplyAD done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaNMCBal::multiplyInverseAD(const Increment & dxin, Increment & dxout) const {
  util::Timer timer(classname(), "multiplyInverseAD");
  oops::Log::trace() << classname() << " multiplyInverseAD starting" << std::endl;
  fv3jedi_linvarcha_nmcbal_multiplyinverseadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                                      dxin.toFortran(), dxout.toFortran());
  oops::Log::trace() << classname() << " multiplyInverseAD starting" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaNMCBal::print(std::ostream & os) const {
  os << classname() << " variable change";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
