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
#include "fv3jedi/State/State.h"
#include "fv3jedi/Tlm/Tlm.h"
#include "fv3jedi/Tlm/Tlm.interface.h"
#include "fv3jedi/Tlm/Traj.interface.h"
#include "fv3jedi/Utilities/Traits.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------
static oops::interface::LinearModelMaker<Traits, Tlm> makerTLM_("FV3JEDITLM");
// -------------------------------------------------------------------------------------------------
Tlm::Tlm(const Geometry & resol, const eckit::Configuration & config)
  : keySelf_(0), tstep_(), trajmap_(), linvars_(), an2model_(), finalVars_()
{
  oops::Log::trace() << "Tlm::Tlm starting" << std::endl;

  // Store time step
  tstep_ = util::Duration(config.getString("tstep"));

  oops::Variables tlvars(config, "tlm variables");
  linvars_ = oops::Variables(resol.fieldsMetaData().getLongNameFromAnyName(tlvars));

  eckit::LocalConfiguration linVarChangeConfig;
  linVarChangeConfig.set("linear variable change name", "Analysis2Model");
  an2model_.reset(new LinearVariableChange(resol, linVarChangeConfig));

  // Implementation
  fv3jedi_tlm_create_f90(keySelf_, resol.toFortran(), config);

  oops::Log::trace() << "Tlm::Tlm done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
Tlm::~Tlm() {
  oops::Log::trace() << "Tlm::~Tlm starting" << std::endl;

  // Implementation
  fv3jedi_tlm_delete_f90(keySelf_);

  // Clear trajecotory
  for (trajIter jtra = trajmap_.begin(); jtra != trajmap_.end(); ++jtra) {
    fv3jedi_traj_wipe_f90(jtra->second);
  }
  oops::Log::trace() << "Tlm::~Tlm done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void Tlm::setTrajectory(const State & xx, State & xlr, const ModelBias & bias) {
  oops::Log::trace() << "Tlm::setTrajectory starting" << std::endl;

  // Interpolate to resolution of the trajectory
  xlr.changeResolution(xx);

  an2model_->changeVarTraj(xlr, linvars_);

  // Set trajectory
  int keyTraj = 0;
  fv3jedi_traj_set_f90(keyTraj, xlr.toFortran());
  ASSERT(keyTraj != 0);
  trajmap_[xx.validTime()] = keyTraj;

  oops::Log::trace() << "Tlm::setTrajectory done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void Tlm::initializeTL(Increment & dx) const {
  oops::Log::trace() << "Tlm::initializeTL starting" << std::endl;

  // Get traj index
  trajICst itra = trajmap_.find(dx.validTime());

  // Check traj is available
  if (itra == trajmap_.end()) {
    oops::Log::error() << "Tlm: trajectory not available at time " << dx.validTime() << std::endl;
    ABORT("Tlm: trajectory not available");
  }

  ASSERT_MSG(!finalVars_, "finalVars_ should always be null when calling initializeTL");
  if (!(linvars_ <= dx.variables())) {
    finalVars_.reset(new oops::Variables(dx.variables()));
    an2model_->changeVarTL(dx, linvars_);
  }

  // Implementation
  fv3jedi_tlm_initialize_tl_f90(keySelf_, dx.toFortran(), itra->second);

  oops::Log::trace() << "Tlm::initializeTL done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void Tlm::stepTL(Increment & dx, const ModelBiasIncrement &) const {
  oops::Log::trace() << "Tlm::stepTL starting" << std::endl;

  // Get traj index
  trajICst itra = trajmap_.find(dx.validTime());

  // Check traj is available
  if (itra == trajmap_.end()) {
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

  if (finalVars_) {
    const bool force_varchange = true;
    an2model_->changeVarInverseTL(dx, *finalVars_, force_varchange);
    finalVars_.reset(nullptr);  // reset to null for next initializeTL
  }

  oops::Log::trace() << "Tlm::finalizeTL done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void Tlm::initializeAD(Increment & dx) const {
  oops::Log::trace() << "Tlm::initializeAD starting" << std::endl;

  // Get traj index
  trajICst itra = trajmap_.find(dx.validTime());

  // Check traj is available
  if (itra == trajmap_.end()) {
    oops::Log::error() << "Tlm: trajectory not available at time " << dx.validTime() << std::endl;
    ABORT("Tlm: trajectory not available");
  }

  ASSERT_MSG(!finalVars_, "finalVars_ should always be null when calling initializeAD");
  if (!(linvars_ <= dx.variables())) {
    finalVars_.reset(new oops::Variables(dx.variables()));
    an2model_->changeVarInverseAD(dx, linvars_);
  }

  // Implementation
  fv3jedi_tlm_initialize_ad_f90(keySelf_, dx.toFortran(), itra->second);

  oops::Log::trace() << "Tlm::initializeAD done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void Tlm::stepAD(Increment & dx, ModelBiasIncrement &) const {
  oops::Log::trace() << "Tlm::stepAD starting" << std::endl;

  // Tick increment clock (backwards)
  dx.validTime() -= tstep_;

  // Get traj index
  trajICst itra = trajmap_.find(dx.validTime());

  // Check traj is available
  if (itra == trajmap_.end()) {
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

  if (finalVars_) {
    const bool force_varchange = true;
    an2model_->changeVarAD(dx, *finalVars_, force_varchange);
    finalVars_.reset(nullptr);  // reset to null for next initializeAD
  }

  oops::Log::trace() << "Tlm::finalizeAD done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void Tlm::print(std::ostream & os) const {
  oops::Log::trace() << "Tlm::print starting" << std::endl;

  // Print information about the Tlm object
  os << "FV3JEDI TLM Trajectory, nstep=" << trajmap_.size() << std::endl;
  typedef std::map< util::DateTime, int >::const_iterator trajICst;
  if (trajmap_.size() > 0) {
    os << "FV3JEDI TLM Trajectory: times are:";
    for (trajICst jtra = trajmap_.begin(); jtra != trajmap_.end(); ++jtra) {
      os << "  " << jtra->first;
    }
  }

  oops::Log::trace() << "Tlm::print done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
