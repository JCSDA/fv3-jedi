/*
 * (C) Copyright 2017-2018  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/VariableChanges/Control2Analysis/LinVarChaC2A.h"

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"
#include "oops/util/Logger.h"
#include "LinVarChaC2A.interface.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
LinVarChaC2A::LinVarChaC2A(const State & bg, const State & fg, const Geometry & resol,
                           const eckit::Configuration & conf):
  geom_(new Geometry(resol))
{
  const eckit::Configuration * configc = &conf;
  fv3jedi_linvarcha_c2a_create_f90(keyFtnConfig_, geom_->toFortran(), bg.toFortran(),
                                   fg.toFortran(), &configc);
  oops::Log::trace() << "LinVarChaC2A created" << std::endl;
}
// -------------------------------------------------------------------------------------------------
LinVarChaC2A::~LinVarChaC2A() {
  fv3jedi_linvarcha_c2a_delete_f90(keyFtnConfig_);
  oops::Log::trace() << "LinVarChaC2A destructed" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void LinVarChaC2A::multiply(const Increment & dxc, Increment & dxa) const {
  oops::Log::trace() << "LinVarChaC2A multiply" << std::endl;
  fv3jedi_linvarcha_c2a_multiply_f90(keyFtnConfig_, geom_->toFortran(),
                                     dxc.toFortran(), dxa.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinVarChaC2A::multiplyInverse(const Increment & dxa, Increment & dxc) const {
  oops::Log::trace() << "LinVarChaC2A multiplyInverse" << std::endl;
  fv3jedi_linvarcha_c2a_multiplyinverse_f90(keyFtnConfig_, geom_->toFortran(),
                                            dxa.toFortran(), dxc.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinVarChaC2A::multiplyAD(const Increment & dxa, Increment & dxc) const {
  oops::Log::trace() << "LinVarChaC2A multiplyAD" << std::endl;
  fv3jedi_linvarcha_c2a_multiplyadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                            dxa.toFortran(), dxc.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinVarChaC2A::multiplyInverseAD(const Increment & dxc, Increment & dxa) const {
  oops::Log::trace() << "LinVarChaC2A multiplyInverseAD" << std::endl;
  fv3jedi_linvarcha_c2a_multiplyinverseadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                                   dxc.toFortran(), dxa.toFortran());
}
// -------------------------------------------------------------------------------------------------
void LinVarChaC2A::print(std::ostream & os) const {
  os << "FV3JEDI change variable";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
