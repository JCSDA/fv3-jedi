/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Tlm/TlmId.h"
#include "fv3jedi/Utilities/Traits.h"
#include "fv3jedi/Utilities/Utilities.h"

namespace fv3jedi {

// -----------------------------------------------------------------------------
static oops::LinearModelMaker<Traits, TlmId>
                                      makerIdTLM_("FV3JEDIIdTLM");
// -----------------------------------------------------------------------------
TlmId::TlmId(const Geometry & resol,
                            const eckit::Configuration & tlConf)
  : keyConfig_(0), tstep_(), resol_(resol), linvars_(tlConf)
{
  tstep_ = util::Duration(tlConf.getString("tstep"));
  oops::Log::trace() << "TlmId created" << std::endl;
}
// -----------------------------------------------------------------------------
TlmId::~TlmId() {
  oops::Log::trace() << "TlmId destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmId::setTrajectory(const State &, State &,
                                 const ModelBias &) {}
// -----------------------------------------------------------------------------
void TlmId::initializeTL(Increment & dx) const {
  oops::Log::debug() << "TlmId::initializTL" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmId::stepTL(Increment & dx,
                             const ModelBiasIncrement &) const {
  dx.updateTime(tstep_);
}
// -----------------------------------------------------------------------------
void TlmId::finalizeTL(Increment & dx) const {
  oops::Log::debug() << "TlmId::finalizeTL" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmId::initializeAD(Increment & dx) const {
  oops::Log::debug() << "TlmId::initializAD" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmId::stepAD(Increment & dx,
                               ModelBiasIncrement &) const {
  dx.updateTime(-tstep_);
}
// -----------------------------------------------------------------------------
void TlmId::finalizeAD(Increment & dx) const {
  oops::Log::debug() << "TlmId::finalizeAD" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmId::print(std::ostream & os) const {
  os << "FV3JEDI IdTLM" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
