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

#include "src/Tlm/TlmIdFV3JEDI.h"
#include "GeometryFV3JEDI.h"
#include "FV3JEDITraits.h"
#include "IncrementFV3JEDI.h"
#include "TlmFV3JEDIFortran.h"
#include "StateFV3JEDI.h"
#include "UtilitiesFV3JEDI.h"

namespace fv3jedi {

// -----------------------------------------------------------------------------
static oops::LinearModelMaker<FV3JEDITraits, TlmIdFV3JEDI>
                                      makerFV3JEDIIdTLM_("FV3JEDIIdTLM");
// -----------------------------------------------------------------------------
TlmIdFV3JEDI::TlmIdFV3JEDI(const GeometryFV3JEDI & resol,
                            const eckit::Configuration & tlConf)
  : keyConfig_(0), tstep_(), resol_(resol),
     linvars_(std::vector<std::string>{"ua", "va", "t", "ps",
                                       "q", "qi", "ql", "o3"})
{
  tstep_ = util::Duration(tlConf.getString("tstep"));

  const eckit::Configuration * configc = &tlConf;
  stageFv3Files(tlConf);
  fv3jedi_tlm_create_f90(&configc, resol_.toFortran(), keyConfig_);
  removeFv3Files();

  oops::Log::trace() << "TlmIdFV3JEDI created" << std::endl;
}
// -----------------------------------------------------------------------------
TlmIdFV3JEDI::~TlmIdFV3JEDI() {
  fv3jedi_tlm_delete_f90(keyConfig_);
  oops::Log::trace() << "TlmIdFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdFV3JEDI::setTrajectory(const StateFV3JEDI &, StateFV3JEDI &,
                                 const ModelBiasFV3JEDI &) {}
// -----------------------------------------------------------------------------
void TlmIdFV3JEDI::initializeTL(IncrementFV3JEDI & dx) const {
  oops::Log::debug() << "TlmIdFV3JEDI::initializTL" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdFV3JEDI::stepTL(IncrementFV3JEDI & dx,
                             const ModelBiasIncrementFV3JEDI &) const {
  dx.updateTime(tstep_);
}
// -----------------------------------------------------------------------------
void TlmIdFV3JEDI::finalizeTL(IncrementFV3JEDI & dx) const {
  oops::Log::debug() << "TlmIdFV3JEDI::finalizeTL" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdFV3JEDI::initializeAD(IncrementFV3JEDI & dx) const {
  oops::Log::debug() << "TlmIdFV3JEDI::initializAD" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdFV3JEDI::stepAD(IncrementFV3JEDI & dx,
                               ModelBiasIncrementFV3JEDI &) const {
  dx.updateTime(-tstep_);
}
// -----------------------------------------------------------------------------
void TlmIdFV3JEDI::finalizeAD(IncrementFV3JEDI & dx) const {
  oops::Log::debug() << "TlmIdFV3JEDI::finalizeAD" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdFV3JEDI::print(std::ostream & os) const {
  os << "FV3JEDI IdTLM" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
