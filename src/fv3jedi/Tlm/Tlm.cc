/*
 * (C) Copyright 2017-2020 UCAR
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

// -------------------------------------------------------------------------------------------------
static oops::LinearModelMaker<Traits, Tlm> makerTLM_("FV3JEDITLM");
// -------------------------------------------------------------------------------------------------
Tlm::Tlm(const Geometry & resol,
                        const eckit::Configuration & tlConf)
  : keySelf_(0), tstep_(), traj_(),
    lrmodel_(resol, eckit::LocalConfiguration(tlConf, "trajectory")),
    linvars_(tlConf, "tlm variables")
{
  oops::Log::trace() << "Tlm::Tlm starting" << std::endl;

  // Store time step
  tstep_ = util::Duration(tlConf.getString("tstep"));

  // Pointer to configuration
  const eckit::Configuration * configc = &tlConf;

  // Implementation
  stageFv3Files(tlConf, resol.getComm());
  fv3jedi_tlm_create_f90(keySelf_, resol.toFortran(), &configc);
  removeFv3Files(resol.getComm());

  oops::Log::trace() << "Tlm::Tlm done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
Tlm::~Tlm() {
  oops::Log::trace() << "Tlm::~Tlm starting" << std::endl;

  // Implementation
  fv3jedi_tlm_delete_f90(keySelf_);

  // Clear trajecotory
  for (trajIter jtra = traj_.begin(); jtra != traj_.end(); ++jtra) {
    fv3jedi_traj_wipe_f90(jtra->second);
  }
  oops::Log::trace() << "Tlm::~Tlm done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void Tlm::setTrajectory(const State & xx, State & xlr, const ModelBias & bias) {
  oops::Log::trace() << "Tlm::setTrajectory starting" << std::endl;

  // Interpolate to resolution of the trajectory
  xlr.changeResolution(xx);

  // Save position in map
  int ftraj = lrmodel_.saveTrajectory(xlr, bias);
  traj_[xx.validTime()] = ftraj;

  // Print method needs implementing
  // std::vector<double> zstat(15);
  // fv3jedi_traj_minmaxrms_f90(ftraj, zstat[0]);
  // oops::Log::debug() << "Tlm trajectory at time " << xx.validTime() << std::endl;
  // for (unsigned int jj = 0; jj < 5; ++jj) {
  //   oops::Log::debug() << "  Min=" << zstat[3*jj] << ", Max=" << zstat[3*jj+1]
  //                      << ", RMS=" << zstat[3*jj+2] << std::endl;
  // }

  oops::Log::trace() << "Tlm::setTrajectory done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void Tlm::initializeTL(Increment & dx) const {
  oops::Log::trace() << "Tlm::initializeTL starting" << std::endl;

  // Implementation
  fv3jedi_tlm_initialize_tl_f90(keySelf_, dx.toFortran());

  oops::Log::trace() << "Tlm::initializeTL done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void Tlm::stepTL(Increment & dx, const ModelBiasIncrement &) const {
  oops::Log::trace() << "Tlm::stepTL starting" << std::endl;

  // Get traj index
  trajICst itra = traj_.find(dx.validTime());

  // Check traj is available
  if (itra == traj_.end()) {
    oops::Log::error() << "Tlm: trajectory not available at time " << dx.validTime() << std::endl;
    ABORT("Tlm: trajectory not available");
  }

  // Implementation
  fv3jedi_tlm_step_tl_f90(keySelf_, dx.toFortran(), itra->second);

  // Tick increment clock
  dx.validTime() += tstep_;

  oops::Log::trace() << "Tlm::stepTL done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void Tlm::finalizeTL(Increment & dx) const {
  oops::Log::trace() << "Tlm::finalizeTL starting" << std::endl;

  // Implementation
  fv3jedi_tlm_finalize_tl_f90(keySelf_, dx.toFortran());

  oops::Log::trace() << "Tlm::finalizeTL done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void Tlm::initializeAD(Increment & dx) const {
  oops::Log::trace() << "Tlm::initializeAD starting" << std::endl;

  // Implementation
  fv3jedi_tlm_initialize_ad_f90(keySelf_, dx.toFortran());

  oops::Log::trace() << "Tlm::initializeAD done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void Tlm::stepAD(Increment & dx, ModelBiasIncrement &) const {
  oops::Log::trace() << "Tlm::stepAD starting" << std::endl;

  // Tick increment clock (backwards)
  dx.validTime() -= tstep_;

  // Get traj index
  trajICst itra = traj_.find(dx.validTime());

  // Check traj is available
  if (itra == traj_.end()) {
    oops::Log::error() << "Tlm: trajectory not available at time " << dx.validTime() << std::endl;
    ABORT("Tlm: trajectory not available");
  }

  // Implementation
  fv3jedi_tlm_step_ad_f90(keySelf_, dx.toFortran(), itra->second);

  oops::Log::trace() << "Tlm::stepAD done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void Tlm::finalizeAD(Increment & dx) const {
  oops::Log::trace() << "Tlm::finalizeAD starting" << std::endl;

  // Implementation
  fv3jedi_tlm_finalize_ad_f90(keySelf_, dx.toFortran());

  oops::Log::trace() << "Tlm::finalizeAD done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void Tlm::print(std::ostream & os) const {
  oops::Log::trace() << "Tlm::print starting" << std::endl;

  // Print information about the Tlm object
  os << "FV3JEDI TLM Trajectory, nstep=" << traj_.size() << std::endl;
  typedef std::map< util::DateTime, int >::const_iterator trajICst;
  if (traj_.size() > 0) {
    os << "FV3JEDI TLM Trajectory: times are:";
    for (trajICst jtra = traj_.begin(); jtra != traj_.end(); ++jtra) {
      os << "  " << jtra->first;
    }
  }

  oops::Log::trace() << "Tlm::print done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
