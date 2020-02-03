/*
 * (C) Copyright 2017-2018  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/VariableChanges/NMCBalance/LinVarChaNMCBal.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"
#include "oops/util/Logger.h"
#include "LinVarChaNMCBal.interface.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
LinVarChaNMCBal::LinVarChaNMCBal(const State & bg, const State & fg, const Geometry & resol,
                                 const eckit::Configuration & conf):
  geom_(new Geometry(resol))
{
  const eckit::Configuration * configc = &conf;
  fv3jedi_linvarcha_nmcbal_create_f90(keyFtnConfig_, geom_->toFortran(), bg.toFortran(),
                                      fg.toFortran(), &configc);
  oops::Log::trace() << "LinVarChaNMCBal created" << std::endl;
}
// -------------------------------------------------------------------------------------------------
LinVarChaNMCBal::~LinVarChaNMCBal() {
  fv3jedi_linvarcha_nmcbal_delete_f90(keyFtnConfig_);
  oops::Log::trace() << "LinVarChaNMCBal destructed" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaNMCBal::multiply(const Increment & dxu, Increment & dxb) const {
  oops::Log::trace() << "LinVarChaNMCBal multiply" << std::endl;
  fv3jedi_linvarcha_nmcbal_multiply_f90(keyFtnConfig_, geom_->toFortran(),
                                        dxu.toFortran(), dxb.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinVarChaNMCBal::multiplyInverse(const Increment & dxb, Increment & dxu) const {
  oops::Log::trace() << "LinVarChaNMCBal multiplyInverse" << std::endl;
  fv3jedi_linvarcha_nmcbal_multiplyinverse_f90(keyFtnConfig_, geom_->toFortran(),
                                               dxb.toFortran(), dxu.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinVarChaNMCBal::multiplyAD(const Increment & dxb, Increment & dxu) const {
  oops::Log::trace() << "LinVarChaNMCBal multiplyAD" << std::endl;
  fv3jedi_linvarcha_nmcbal_multiplyadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                               dxb.toFortran(), dxu.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinVarChaNMCBal::multiplyInverseAD(const Increment & dxu, Increment & dxb) const {
  oops::Log::trace() << "LinVarChaNMCBal multiplyInverseAD" << std::endl;
  fv3jedi_linvarcha_nmcbal_multiplyinverseadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                                      dxu.toFortran(), dxb.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinVarChaNMCBal::print(std::ostream & os) const {
  os << "FV3JEDI change variable";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
