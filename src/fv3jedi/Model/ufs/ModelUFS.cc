/*
 * (C) Copyright 2020 NOAA
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "ModelUFS.interface.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Model/ufs/ModelUFS.h"
#include "fv3jedi/ModelBias/ModelBias.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Utilities/Utilities.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static oops::ModelMaker<Traits, ModelUFS> makermodel_("UFS");
// -------------------------------------------------------------------------------------------------
ModelUFS::ModelUFS(const Geometry & resol, const eckit::Configuration & mconf)
  : keyConfig_(0), tstep_(0), geom_(resol), vars_(mconf, "model variables") {
  oops::Log::trace() << "ModelUFS::ModelUFS starting" << std::endl;
  char tmpdir[10000];
  getcwd(tmpdir, 10000);
  tstep_ = util::Duration(mconf.getString("tstep"));
  strcpy(ufsdir_, mconf.getString("ufs_run_directory").c_str());
  const eckit::Configuration * configc = &mconf;
  chdir(ufsdir_);
  fv3jedi_ufs_create_f90(keyConfig_, &configc, geom_.toFortran());
  oops::Log::trace() << "ModelUFS::ModelUFS done" << std::endl;
  chdir(tmpdir);
}
// -------------------------------------------------------------------------------------------------
ModelUFS::~ModelUFS() {
  oops::Log::trace() << "ModelUFS::~ModelUFS starting" << std::endl;
  fv3jedi_ufs_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelUFS::~ModelUFS done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelUFS::initialize(State & xx) const {
  oops::Log::trace() << "ModelUFS::initialize starting" << std::endl;
  oops::Log::trace() << "ModelUFS::cd to " << ufsdir_ << std::endl;

  chdir(ufsdir_);
  util::DateTime * dtp = &xx.validTime();
  fv3jedi_ufs_initialize_f90(keyConfig_, xx.toFortran(), &dtp);
  oops::Log::trace() << "ModelUFS::initialize done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelUFS::step(State & xx, const ModelBias &) const
{
  oops::Log::trace() << "ModelUFS::step starting" << std::endl;
  oops::Log::trace() << "ModelUFS::cd to " << ufsdir_ << std::endl;
  char tmpdir[10000];
  chdir(ufsdir_);
  getcwd(tmpdir, 10000);
  oops::Log::trace() << "ModelUFS::PWD is  " << tmpdir << std::endl;

  util::DateTime start = xx.validTime();
  util::DateTime * dtp1 = &start;
  oops::Log::trace() << "Model start time is " << xx.validTime() << std::endl;
  oops::Log::trace() << "Forecast time step is " << tstep_ << std::endl;
  xx.validTime() += tstep_;
  util::DateTime * dtp2 = &xx.validTime();
  fv3jedi_ufs_step_f90(keyConfig_, xx.toFortran(), &dtp1, &dtp2);
  oops::Log::trace() << "ModelUFS::step done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelUFS::finalize(State & xx) const {
  oops::Log::trace() << "ModelUFS::finalize starting" << std::endl;
  util::DateTime * dtp = &xx.validTime();
  fv3jedi_ufs_finalize_f90(keyConfig_, xx.toFortran(), &dtp);
  oops::Log::trace() << "ModelUFS::finalize done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
int ModelUFS::saveTrajectory(State & xx, const ModelBias &) const {
  ABORT("Model:UFS should not be used for the trajecotry");
}
// -------------------------------------------------------------------------------------------------
void ModelUFS::print(std::ostream & os) const {
  os << "ModelUFS::print not implemented";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
