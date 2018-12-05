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

#include "ModelNUOPCFV3JEDIFortran.h"
#include "GeometryFV3JEDI.h"
#include "ModelBiasFV3JEDI.h"
#include "src/Model/nuopc/ModelNUOPCFV3JEDI.h"
#include "StateFV3JEDI.h"
#include "UtilitiesFV3JEDI.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
static oops::ModelMaker<FV3JEDITraits, ModelNUOPCFV3JEDI> makermodel_("NUOPC");
// -----------------------------------------------------------------------------
ModelNUOPCFV3JEDI::ModelNUOPCFV3JEDI(const GeometryFV3JEDI & resol,
                            const eckit::Configuration & model)
  : keyConfig_(0), tstep_(0), geom_(resol),
  vars_(std::vector<std::string>{"ud", "vd", "ua", "va", "t", "delp",
                                 "q", "qi", "ql", "o3"})
{
  oops::Log::trace() << "ModelNUOPCFV3JEDI::ModelNUOPCFV3JEDI" << std::endl;
  tstep_ = util::Duration(model.getString("tstep"));
  const eckit::Configuration * configc = &model;
  stageFv3Files(model);
  fv3jedi_nuopc_create_f90(&configc, geom_.toFortran(), keyConfig_);
  removeFv3Files();
  oops::Log::trace() << "ModelNUOPCFV3JEDI created" << std::endl;
}
// -----------------------------------------------------------------------------
ModelNUOPCFV3JEDI::~ModelNUOPCFV3JEDI() {
  fv3jedi_nuopc_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelNUOPCFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelNUOPCFV3JEDI::initialize(StateFV3JEDI & xx) const {
  util::DateTime * dtp = &xx.validTime();
  fv3jedi_nuopc_initialize_f90(keyConfig_, xx.toFortran(), &dtp);
  oops::Log::debug() << "ModelNUOPCFV3JEDI::initialize" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelNUOPCFV3JEDI::step(StateFV3JEDI & xx, const ModelBiasFV3JEDI &) const {
  util::DateTime start = xx.validTime();
  util::DateTime * dtp1 = &start;
  xx.validTime() += tstep_;
  util::DateTime * dtp2 = &xx.validTime();
  fv3jedi_nuopc_step_f90(keyConfig_, xx.toFortran(), &dtp1, &dtp2);
  oops::Log::debug() << "ModelNUOPCFV3JEDI::step" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelNUOPCFV3JEDI::finalize(StateFV3JEDI & xx) const {
  util::DateTime * dtp = &xx.validTime();
  fv3jedi_nuopc_finalize_f90(keyConfig_, xx.toFortran(), &dtp);
  oops::Log::debug() << "ModelNUOPCFV3JEDI::finalize" << std::endl;
}
// -----------------------------------------------------------------------------
int ModelNUOPCFV3JEDI::saveTrajectory(StateFV3JEDI & xx,
                                 const ModelBiasFV3JEDI &) const {
  int ftraj = 0;
  fv3jedi_traj_prop_f90(keyConfig_, xx.toFortran(), ftraj);
  ASSERT(ftraj != 0);
  return ftraj;
}
// -----------------------------------------------------------------------------
void ModelNUOPCFV3JEDI::print(std::ostream & os) const {
  os << "ModelNUOPCFV3JEDI::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
