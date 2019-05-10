/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "ModelPseudoFV3JEDIFortran.h"
#include "GeometryFV3JEDI.h"
#include "ModelBiasFV3JEDI.h"
#include "src/Model/pseudo/ModelPseudoFV3JEDI.h"
#include "StateFV3JEDI.h"
#include "UtilitiesFV3JEDI.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
static oops::ModelMaker<FV3JEDITraits, ModelPseudoFV3JEDI>
                        makermodel_("PSEUDO");
// -----------------------------------------------------------------------------
ModelPseudoFV3JEDI::ModelPseudoFV3JEDI(const GeometryFV3JEDI & resol,
                            const eckit::Configuration & mconf)
  : keyConfig_(0), tstep_(0), geom_(resol), vars_(mconf)
{
  oops::Log::trace() << "ModelPseudoFV3JEDI::ModelPseudoFV3JEDI" << std::endl;
  tstep_ = util::Duration(mconf.getString("tstep"));
  const eckit::Configuration * configc = &mconf;
  fv3jedi_pseudo_create_f90(&configc, geom_.toFortran(), keyConfig_);
  if (mconf.has("RunStageCheck")) {
    runstagecheck_ = mconf.getInt("RunStageCheck");
  }
  oops::Log::trace() << "ModelPseudoFV3JEDI created" << std::endl;
}
// -----------------------------------------------------------------------------
ModelPseudoFV3JEDI::~ModelPseudoFV3JEDI() {
  fv3jedi_pseudo_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelPseudoFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelPseudoFV3JEDI::initialize(StateFV3JEDI & xx) const {
  if (runstage_) {
    fv3jedi_pseudo_initialize_f90(keyConfig_, xx.toFortran());
  }
  oops::Log::debug() << "ModelPseudoFV3JEDI::initialize" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelPseudoFV3JEDI::step(StateFV3JEDI & xx,
                              const ModelBiasFV3JEDI &) const {
  xx.validTime() += tstep_;
  util::DateTime * dtp = &xx.validTime();
  if (runstage_) {
    fv3jedi_pseudo_step_f90(keyConfig_, xx.toFortran(),
                            geom_.toFortran(), &dtp);
  } else {
    oops::Log::warning() << "Pseudo model has already run through once so not"
                            "re-reading, just ticking the clock." << std::endl;
  }
  oops::Log::debug() << "ModelPseudoFV3JEDI::step" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelPseudoFV3JEDI::finalize(StateFV3JEDI & xx) const {
  if (runstage_) {
    fv3jedi_pseudo_finalize_f90(keyConfig_, xx.toFortran());
  }
  if (runstagecheck_ == 1) {runstage_ = false;}
  oops::Log::debug() << "ModelPseudoFV3JEDI::finalize" << std::endl;
}
// -----------------------------------------------------------------------------
int ModelPseudoFV3JEDI::saveTrajectory(StateFV3JEDI & xx,
                                 const ModelBiasFV3JEDI &) const {
  int ftraj = 0;
  fv3jedi_traj_prop_f90(keyConfig_, xx.toFortran(), ftraj);
  ASSERT(ftraj != 0);
  return ftraj;
}
// -----------------------------------------------------------------------------
void ModelPseudoFV3JEDI::print(std::ostream & os) const {
  os << "ModelPseudoFV3JEDI::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
