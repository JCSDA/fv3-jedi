/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <vector>

#include "TlmFV3JEDI.h"

#include "eckit/config/LocalConfiguration.h"
#include "Fortran.h"
#include "GeometryFV3JEDI.h"
#include "IncrementFV3JEDI.h"
#include "ModelFV3JEDI.h"
#include "StateFV3JEDI.h"
#include "FV3JEDITraits.h"
#include "util/DateTime.h"
#include "util/Logger.h"
#include "util/abor1_cpp.h"
#include "UtilitiesFV3JEDI.h"

namespace fv3jedi {

// -----------------------------------------------------------------------------
static oops::LinearModelMaker<FV3JEDITraits, TlmFV3JEDI> makerFV3JEDITLM_("FV3JEDITLM");
// -----------------------------------------------------------------------------
TlmFV3JEDI::TlmFV3JEDI(const GeometryFV3JEDI & resol, const eckit::Configuration & tlConf)
  : keyConfig_(0), tstep_(), resol_(resol), traj_(),
    lrmodel_(resol_, eckit::LocalConfiguration(tlConf, "trajectory"))
{
  tstep_ = util::Duration(tlConf.getString("tstep"));

  const eckit::Configuration * configc = &tlConf;
  stageFv3FilesPert(tlConf);
  fv3jedi_model_setup_f90(&configc, resol_.toFortran(), keyConfig_);
  removeFv3FilesPert();
  oops::Log::trace() << "TlmFV3JEDI created" << std::endl;
}
// -----------------------------------------------------------------------------
TlmFV3JEDI::~TlmFV3JEDI() {
  fv3jedi_model_delete_f90(keyConfig_);
  for (trajIter jtra = traj_.begin(); jtra != traj_.end(); ++jtra) {
    fv3jedi_model_wipe_traj_f90(jtra->second);
  }
  oops::Log::trace() << "TlmFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmFV3JEDI::setTrajectory(const StateFV3JEDI & xx, StateFV3JEDI & xlr, const ModelBiasFV3JEDI & bias) {
// StateFV3JEDI xlr(resol_, xx);
  xlr.changeResolution(xx);
  int ftraj = lrmodel_.saveTrajectory(xlr, bias);
  traj_[xx.validTime()] = ftraj;

// should be in print method
  std::vector<double> zstat(15);
//  fv3jedi_traj_minmaxrms_f90(ftraj, zstat[0]);
  oops::Log::debug() << "TlmFV3JEDI trajectory at time " << xx.validTime() << std::endl;
  for (unsigned int jj = 0; jj < 5; ++jj) {
    oops::Log::debug() << "  Min=" << zstat[3*jj] << ", Max=" << zstat[3*jj+1]
                       << ", RMS=" << zstat[3*jj+2] << std::endl;
  }
// should be in print method
}
// -----------------------------------------------------------------------------
void TlmFV3JEDI::initializeTL(IncrementFV3JEDI & dx) const {
  fv3jedi_model_prepare_integration_tl_f90(keyConfig_, dx.fields().toFortran());
  oops::Log::debug() << "TlmFV3JEDI::initializeTL" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmFV3JEDI::stepTL(IncrementFV3JEDI & dx, const ModelBiasIncrementFV3JEDI &) const {
  trajICst itra = traj_.find(dx.validTime());
  if (itra == traj_.end()) {
    oops::Log::error() << "TlmFV3JEDI: trajectory not available at time " << dx.validTime() << std::endl;
    ABORT("TlmFV3JEDI: trajectory not available");
  }
  oops::Log::debug() << "TlmFV3JEDI::stepTL fields in" << dx.fields() << std::endl;
  fv3jedi_model_propagate_tl_f90(keyConfig_, dx.fields().toFortran(), itra->second);
  oops::Log::debug() << "TlmFV3JEDI::stepTL fields out" << dx.fields() << std::endl;
  dx.validTime() += tstep_;
}
// -----------------------------------------------------------------------------
void TlmFV3JEDI::finalizeTL(IncrementFV3JEDI & dx) const {
  oops::Log::debug() << "TlmFV3JEDI::finalizeTL" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmFV3JEDI::initializeAD(IncrementFV3JEDI & dx) const {
  oops::Log::debug() << "TlmFV3JEDI::initializeAD" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmFV3JEDI::stepAD(IncrementFV3JEDI & dx, ModelBiasIncrementFV3JEDI &) const {
  dx.validTime() -= tstep_;
  trajICst itra = traj_.find(dx.validTime());
  if (itra == traj_.end()) {
    oops::Log::error() << "TlmFV3JEDI: trajectory not available at time " << dx.validTime() << std::endl;
    ABORT("TlmFV3JEDI: trajectory not available");
  }
  oops::Log::debug() << "TlmFV3JEDI::stepAD fields in" << dx.fields() << std::endl;
  fv3jedi_model_propagate_ad_f90(keyConfig_, dx.fields().toFortran(), itra->second);
  oops::Log::debug() << "TlmFV3JEDI::stepAD fields out" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmFV3JEDI::finalizeAD(IncrementFV3JEDI & dx) const {
  fv3jedi_model_prepare_integration_ad_f90(keyConfig_, dx.fields().toFortran());
  oops::Log::debug() << "TlmFV3JEDI::finalizeAD" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmFV3JEDI::print(std::ostream & os) const {
  os << "FV3JEDI TLM Trajectory, nstep=" << traj_.size() << std::endl;
  typedef std::map< util::DateTime, int >::const_iterator trajICst;
  if (traj_.size() > 0) {
    os << "FV3JEDI TLM Trajectory: times are:";
    for (trajICst jtra = traj_.begin(); jtra != traj_.end(); ++jtra) {
      os << "  " << jtra->first;
    }
  }
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
