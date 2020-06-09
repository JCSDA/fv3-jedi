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
#include "fv3jedi/Model/traj/ModelTraj.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Tlm/Tlm.h"
#include "fv3jedi/Utilities/Traits.h"
#include "fv3jedi/Utilities/Utilities.h"

namespace fv3jedi {

// -----------------------------------------------------------------------------
static oops::LinearModelMaker<Traits, Tlm>
                                     makerTLM_("FV3JEDITLM");
// -----------------------------------------------------------------------------
Tlm::Tlm(const Geometry & resol,
                        const eckit::Configuration & tlConf)
  : keyConfig_(0), tstep_(), resol_(resol), traj_(),
    lrmodel_(resol_, eckit::LocalConfiguration(tlConf, "trajectory")),
    linvars_(tlConf)
{
  tstep_ = util::Duration(tlConf.getString("tstep"));

  const eckit::Configuration * configc = &tlConf;

  stageFv3Files(tlConf, resol_.getComm());
  fv3jedi_tlm_create_f90(&configc, resol_.toFortran(), keyConfig_, linvars_);
  removeFv3Files(resol_.getComm());
  oops::Log::trace() << "Tlm created" << std::endl;
}
// -----------------------------------------------------------------------------
Tlm::~Tlm() {
  fv3jedi_tlm_delete_f90(keyConfig_);
  for (trajIter jtra = traj_.begin(); jtra != traj_.end(); ++jtra) {
    fv3jedi_traj_wipe_f90(jtra->second);
  }
  oops::Log::trace() << "Tlm destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void Tlm::setTrajectory(const State & xx, State & xlr,
                                const ModelBias & bias) {
// State xlr(resol_, xx);
  xlr.changeResolution(xx);
  int ftraj = lrmodel_.saveTrajectory(xlr, bias);
  traj_[xx.validTime()] = ftraj;

// should be in print method
  std::vector<double> zstat(15);
//  fv3jedi_traj_minmaxrms_f90(ftraj, zstat[0]);
  oops::Log::debug() << "Tlm trajectory at time "
                     << xx.validTime() << std::endl;
  for (unsigned int jj = 0; jj < 5; ++jj) {
    oops::Log::debug() << "  Min=" << zstat[3*jj] << ", Max=" << zstat[3*jj+1]
                       << ", RMS=" << zstat[3*jj+2] << std::endl;
  }
// should be in print method
}
// -----------------------------------------------------------------------------
void Tlm::initializeTL(Increment & dx) const {
  fv3jedi_tlm_initialize_tl_f90(resol_.toFortran(), keyConfig_, dx.toFortran());
  oops::Log::debug() << "Tlm::initializeTL" << std::endl;
}
// -----------------------------------------------------------------------------
void Tlm::stepTL(Increment & dx,
                         const ModelBiasIncrement &) const {
  trajICst itra = traj_.find(dx.validTime());
  if (itra == traj_.end()) {
    oops::Log::error() << "Tlm: trajectory not available at time "
                       << dx.validTime() << std::endl;
    ABORT("Tlm: trajectory not available");
  }
  fv3jedi_tlm_step_tl_f90(resol_.toFortran(), keyConfig_,
                                 dx.toFortran(),
                                  itra->second);
  dx.validTime() += tstep_;
}
// -----------------------------------------------------------------------------
void Tlm::finalizeTL(Increment & dx) const {
  fv3jedi_tlm_finalize_tl_f90(resol_.toFortran(), keyConfig_, dx.toFortran());
  oops::Log::debug() << "Tlm::finalizeTL" << std::endl;
}
// -----------------------------------------------------------------------------
void Tlm::initializeAD(Increment & dx) const {
  fv3jedi_tlm_initialize_ad_f90(resol_.toFortran(), keyConfig_, dx.toFortran());
  oops::Log::debug() << "Tlm::initializeAD" << std::endl;
}
// -----------------------------------------------------------------------------
void Tlm::stepAD(Increment & dx, ModelBiasIncrement &)
    const {
  dx.validTime() -= tstep_;
  trajICst itra = traj_.find(dx.validTime());
  if (itra == traj_.end()) {
    oops::Log::error() << "Tlm: trajectory not available at time "
                       << dx.validTime() << std::endl;
    ABORT("Tlm: trajectory not available");
  }
  fv3jedi_tlm_step_ad_f90(resol_.toFortran(), keyConfig_,
                                 dx.toFortran(),
                                  itra->second);
}
// -----------------------------------------------------------------------------
void Tlm::finalizeAD(Increment & dx) const {
  fv3jedi_tlm_finalize_ad_f90(resol_.toFortran(), keyConfig_, dx.toFortran());
  oops::Log::debug() << "Tlm::finalizeAD" << std::endl;
}
// -----------------------------------------------------------------------------
void Tlm::print(std::ostream & os) const {
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
