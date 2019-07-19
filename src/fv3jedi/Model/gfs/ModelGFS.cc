/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Model/gfs/ModelGFS.h"
#include "fv3jedi/ModelBias/ModelBias.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Utilities/Utilities.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
static oops::ModelMaker<Traits, ModelGFS> makermodel_("GFS");
// -----------------------------------------------------------------------------
ModelGFS::ModelGFS(const Geometry & resol,
                            const eckit::Configuration & mconf)
  : keyConfig_(0), tstep_(0), geom_(resol), vars_(mconf)
{
  oops::Log::trace() << "ModelGFS::ModelGFS" << std::endl;
  tstep_ = util::Duration(mconf.getString("tstep"));
  const eckit::Configuration * configc = &mconf;
  fv3jedi_gfs_create_f90(&configc, geom_.toFortran(), keyConfig_);
  oops::Log::trace() << "ModelGFS created" << std::endl;
}
// -----------------------------------------------------------------------------
ModelGFS::~ModelGFS() {
  fv3jedi_gfs_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelGFS destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelGFS::initialize(State & xx) const {
  fv3jedi_gfs_initialize_f90(keyConfig_, xx.toFortran());
  oops::Log::debug() << "ModelGFS::initialize" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelGFS::step(State & xx, const ModelBias &) const {
  xx.validTime() += tstep_;
  util::DateTime * dtp = &xx.validTime();
  fv3jedi_gfs_step_f90(keyConfig_, xx.toFortran(), &dtp);
  oops::Log::debug() << "ModelGFS::step" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelGFS::finalize(State & xx) const {
  fv3jedi_gfs_finalize_f90(keyConfig_, xx.toFortran());
  oops::Log::debug() << "ModelGFS::finalize" << std::endl;
}
// -----------------------------------------------------------------------------
int ModelGFS::saveTrajectory(State & xx,
                                 const ModelBias &) const {
  ABORT("Model:GFS should not be used for the trajecotry");
}
// -----------------------------------------------------------------------------
void ModelGFS::print(std::ostream & os) const {
  os << "ModelGFS::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
