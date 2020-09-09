/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Model/pseudo/ModelPseudo.h"
#include "fv3jedi/ModelBias/ModelBias.h"
#include "fv3jedi/State/State.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static oops::ModelMaker<Traits, ModelPseudo> makermodel_("PSEUDO");
// -------------------------------------------------------------------------------------------------
ModelPseudo::ModelPseudo(const Geometry & resol, const eckit::Configuration & mconf)
  : keyConfig_(0), tstep_(0), geom_(resol), vars_(mconf, "model variables")
{
  oops::Log::trace() << "ModelPseudo::ModelPseudo" << std::endl;
  tstep_ = util::Duration(mconf.getString("tstep"));
  const eckit::Configuration * configc = &mconf;
  fv3jedi_pseudo_create_f90(&configc, geom_.toFortran(), keyConfig_);
  if (mconf.has("run stage check")) {
    runstagecheck_ = mconf.getInt("run stage check");
  }
  oops::Log::trace() << "ModelPseudo created" << std::endl;
}
// -------------------------------------------------------------------------------------------------
ModelPseudo::~ModelPseudo() {
  fv3jedi_pseudo_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelPseudo destructed" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelPseudo::initialize(State & xx) const {
  if (runstage_) {
    fv3jedi_pseudo_initialize_f90(keyConfig_, xx.toFortran());
  }
  oops::Log::trace() << "ModelPseudo::initialize" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelPseudo::step(State & xx, const ModelBias &) const {
  xx.validTime() += tstep_;
  util::DateTime * dtp = &xx.validTime();
  if (runstage_) {
    fv3jedi_pseudo_step_f90(keyConfig_, xx.toFortran(), geom_.toFortran(), &dtp);
  } else {
    int world_rank = oops::mpi::world().rank();
    if (world_rank == 0) {
      oops::Log::warning() << "Pseudo model has already run through "
                            "once so not re-reading states, just ticking the"
                            " clock." << std::endl;
    }
  }
  oops::Log::trace() << "ModelPseudo::step" << xx.validTime() << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelPseudo::finalize(State & xx) const {
  if (runstage_) {
    fv3jedi_pseudo_finalize_f90(keyConfig_, xx.toFortran());
  }
  if (runstagecheck_ == 1) {runstage_ = false;}
  oops::Log::trace() << "ModelPseudo::finalize" << std::endl;
}
// -------------------------------------------------------------------------------------------------
int ModelPseudo::saveTrajectory(State & xx,
                                 const ModelBias &) const {
  ABORT("Model: pseudo should not be used for the trajecotry");
  return -1;
}
// -------------------------------------------------------------------------------------------------
void ModelPseudo::print(std::ostream & os) const {
  os << "ModelPseudo::print not implemented";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
