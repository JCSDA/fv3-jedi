/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "ModelGFSFV3JEDIFortran.h"
#include "GeometryFV3JEDI.h"
#include "ModelBiasFV3JEDI.h"
#include "src/Model/gfs/ModelGFSFV3JEDI.h"
#include "StateFV3JEDI.h"
#include "UtilitiesFV3JEDI.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
static oops::ModelMaker<FV3JEDITraits, ModelGFSFV3JEDI> makermodel_("GFS");
// -----------------------------------------------------------------------------
ModelGFSFV3JEDI::ModelGFSFV3JEDI(const GeometryFV3JEDI & resol,
                            const eckit::Configuration & model)
  : keyConfig_(0), tstep_(0), geom_(resol),
  vars_(std::vector<std::string>{"ud", "vd", "ua", "va", "t", "delp",
                                 "q", "qi", "ql", "o3"})
{
  oops::Log::trace() << "ModelGFSFV3JEDI::ModelGFSFV3JEDI" << std::endl;
  tstep_ = util::Duration(model.getString("tstep"));
  const eckit::Configuration * configc = &model;
  fv3jedi_gfs_create_f90(&configc, geom_.toFortran(), keyConfig_);
  oops::Log::trace() << "ModelGFSFV3JEDI created" << std::endl;
}
// -----------------------------------------------------------------------------
ModelGFSFV3JEDI::~ModelGFSFV3JEDI() {
  fv3jedi_gfs_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelGFSFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelGFSFV3JEDI::initialize(StateFV3JEDI & xx) const {
  fv3jedi_gfs_initialize_f90(keyConfig_, xx.toFortran());
  oops::Log::debug() << "ModelGFSFV3JEDI::initialize" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelGFSFV3JEDI::step(StateFV3JEDI & xx, const ModelBiasFV3JEDI &) const {
  xx.validTime() += tstep_;
  util::DateTime * dtp = &xx.validTime();
  fv3jedi_gfs_step_f90(keyConfig_, xx.toFortran(), &dtp);
  oops::Log::debug() << "ModelGFSFV3JEDI::step" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelGFSFV3JEDI::finalize(StateFV3JEDI & xx) const {
  fv3jedi_gfs_finalize_f90(keyConfig_, xx.toFortran());
  oops::Log::debug() << "ModelGFSFV3JEDI::finalize" << std::endl;
}
// -----------------------------------------------------------------------------
int ModelGFSFV3JEDI::saveTrajectory(StateFV3JEDI & xx,
                                 const ModelBiasFV3JEDI &) const {
  int ftraj = 0;
  fv3jedi_traj_prop_f90(keyConfig_, xx.toFortran(), ftraj);
  ASSERT(ftraj != 0);
  return ftraj;
}
// -----------------------------------------------------------------------------
void ModelGFSFV3JEDI::print(std::ostream & os) const {
  os << "ModelGFSFV3JEDI::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
