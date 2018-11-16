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

#include "ModelGEOSFV3JEDIFortran.h"
#include "GeometryFV3JEDI.h"
#include "ModelBiasFV3JEDI.h"
#include "src/Model/geos/ModelGEOSFV3JEDI.h"
#include "StateFV3JEDI.h"
#include "UtilitiesFV3JEDI.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
static oops::ModelMaker<FV3JEDITraits, ModelGEOSFV3JEDI> makermodel_("GEOS");
// -----------------------------------------------------------------------------
ModelGEOSFV3JEDI::ModelGEOSFV3JEDI(const GeometryFV3JEDI & resol,
                            const eckit::Configuration & model)
  : keyConfig_(0), tstep_(0), geom_(resol),
  vars_(std::vector<std::string>{"ud", "vd", "ua", "va", "t", "delp",
                                 "q", "qi", "ql", "o3"})
{
  oops::Log::trace() << "ModelGEOSFV3JEDI::ModelGEOSFV3JEDI" << std::endl;
  tstep_ = util::Duration(model.getString("tstep"));
  const eckit::Configuration * configc = &model;
  stageFv3Files(model);
  fv3jedi_geos_create_f90(&configc, geom_.toFortran(), keyConfig_);
  removeFv3Files();
  oops::Log::trace() << "ModelGEOSFV3JEDI created" << std::endl;
}
// -----------------------------------------------------------------------------
ModelGEOSFV3JEDI::~ModelGEOSFV3JEDI() {
  fv3jedi_geos_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelGEOSFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelGEOSFV3JEDI::initialize(StateFV3JEDI & xx) const {
  fv3jedi_geos_initialize_f90(keyConfig_, xx.toFortran());
  oops::Log::debug() << "ModelGEOSFV3JEDI::initialize" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelGEOSFV3JEDI::step(StateFV3JEDI & xx, const ModelBiasFV3JEDI &) const {
  xx.validTime() += tstep_;
  util::DateTime * dtp = &xx.validTime();
  fv3jedi_geos_step_f90(keyConfig_, xx.toFortran(), &dtp);
  oops::Log::debug() << "ModelGEOSFV3JEDI::step" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelGEOSFV3JEDI::finalize(StateFV3JEDI & xx) const {
  fv3jedi_geos_finalize_f90(keyConfig_, xx.toFortran());
  oops::Log::debug() << "ModelGEOSFV3JEDI::finalize" << std::endl;
}
// -----------------------------------------------------------------------------
int ModelGEOSFV3JEDI::saveTrajectory(StateFV3JEDI & xx,
                                 const ModelBiasFV3JEDI &) const {
  int ftraj = 0;
  fv3jedi_traj_prop_f90(keyConfig_, xx.toFortran(), ftraj);
  ASSERT(ftraj != 0);
  return ftraj;
}
// -----------------------------------------------------------------------------
void ModelGEOSFV3JEDI::print(std::ostream & os) const {
  os << "ModelGEOSFV3JEDI::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
