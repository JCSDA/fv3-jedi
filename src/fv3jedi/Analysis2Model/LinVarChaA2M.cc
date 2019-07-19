/*
 * (C) Copyright 2017-2018  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Analysis2Model/LinVarChaA2M.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"
#include "oops/util/Logger.h"
#include "LinVarChaA2M.interface.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
LinVarChaA2M::LinVarChaA2M(const State & bg,
                                   const State & fg,
                                   const Geometry & resol,
                                   const eckit::Configuration & conf):
  geom_(new Geometry(resol))
{
  const eckit::Configuration * configc = &conf;
  fv3jedi_linvarcha_a2m_create_f90(keyFtnConfig_, geom_->toFortran(),
                                   bg.toFortran(),
                                   fg.toFortran(),
                                   &configc);
  oops::Log::trace() << "LinVarChaA2M created" << std::endl;
}
// -----------------------------------------------------------------------------
LinVarChaA2M::~LinVarChaA2M() {
  fv3jedi_linvarcha_a2m_delete_f90(keyFtnConfig_);
  oops::Log::trace() << "ChangeFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void LinVarChaA2M::multiply(const Increment & dxa,
                                         Increment & dxm) const {
  fv3jedi_linvarcha_a2m_multiply_f90(keyFtnConfig_, geom_->toFortran(),
                                     dxa.toFortran(), dxm.toFortran());
  dxm.validTime() = dxa.validTime();
}
// -----------------------------------------------------------------------------
void LinVarChaA2M::multiplyInverse(const Increment & dxm,
                                                Increment & dxa) const {
  fv3jedi_linvarcha_a2m_multiplyinverse_f90(keyFtnConfig_, geom_->toFortran(),
                                            dxm.toFortran(), dxa.toFortran());
  dxa.validTime() = dxm.validTime();
}
// -----------------------------------------------------------------------------
void LinVarChaA2M::multiplyAD(const Increment & dxm,
                                           Increment & dxa) const {
  fv3jedi_linvarcha_a2m_multiplyadjoint_f90(keyFtnConfig_, geom_->toFortran(),
                                            dxm.toFortran(), dxa.toFortran());
  dxa.validTime() = dxm.validTime();
}
// -----------------------------------------------------------------------------
void LinVarChaA2M::multiplyInverseAD(const Increment & dxa,
                                                 Increment & dxm) const {
  fv3jedi_linvarcha_a2m_multiplyinverseadjoint_f90(keyFtnConfig_,
                                                   geom_->toFortran(),
                                                   dxa.toFortran(),
                                                   dxm.toFortran());
  dxm.validTime() = dxa.validTime();
}
// -----------------------------------------------------------------------------
void LinVarChaA2M::print(std::ostream & os) const {
  os << "FV3JEDI change variable";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
