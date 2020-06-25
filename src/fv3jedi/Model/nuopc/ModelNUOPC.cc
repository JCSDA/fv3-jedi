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

#include "ModelNUOPC.interface.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Model/nuopc/ModelNUOPC.h"
#include "fv3jedi/ModelBias/ModelBias.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Utilities/Utilities.h"

namespace fv3jedi {
// -----------------------------------------------------------------------------
static oops::ModelMaker<Traits, ModelNUOPC> makermodel_("NUOPC");
// -----------------------------------------------------------------------------
ModelNUOPC::ModelNUOPC(const Geometry & resol,
                            const eckit::Configuration & model)
  : keyConfig_(0), tstep_(0), geom_(resol),
  vars_(std::vector<std::string>{"ud", "vd", "ua", "va", "t", "delp",
                                 "q", "qi", "ql", "o3mr"})
{
  oops::Log::trace() << "ModelNUOPC::ModelNUOPC" << std::endl;
  tstep_ = util::Duration(model.getString("tstep"));
  const eckit::Configuration * configc = &model;
  fv3jedi_nuopc_create_f90(&configc, geom_.toFortran(), keyConfig_);
  oops::Log::trace() << "ModelNUOPC created" << std::endl;
}
// -----------------------------------------------------------------------------
ModelNUOPC::~ModelNUOPC() {
  fv3jedi_nuopc_delete_f90(keyConfig_);
  oops::Log::trace() << "ModelNUOPC destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelNUOPC::initialize(State & xx) const {
  util::DateTime * dtp = &xx.validTime();
  fv3jedi_nuopc_initialize_f90(keyConfig_, xx.toFortran(), &dtp);
  oops::Log::debug() << "ModelNUOPC::initialize" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelNUOPC::step(State & xx, const ModelBias &) const
{
  util::DateTime start = xx.validTime();
  util::DateTime * dtp1 = &start;
  xx.validTime() += tstep_;
  util::DateTime * dtp2 = &xx.validTime();
  fv3jedi_nuopc_step_f90(keyConfig_, xx.toFortran(), &dtp1, &dtp2);
  oops::Log::debug() << "ModelNUOPC::step" << std::endl;
}
// -----------------------------------------------------------------------------
void ModelNUOPC::finalize(State & xx) const {
  util::DateTime * dtp = &xx.validTime();
  fv3jedi_nuopc_finalize_f90(keyConfig_, xx.toFortran(), &dtp);
  oops::Log::debug() << "ModelNUOPC::finalize" << std::endl;
}
// -----------------------------------------------------------------------------
int ModelNUOPC::saveTrajectory(State & xx,
                                 const ModelBias &) const {
}
// -----------------------------------------------------------------------------
void ModelNUOPC::print(std::ostream & os) const {
  os << "ModelNUOPC::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
