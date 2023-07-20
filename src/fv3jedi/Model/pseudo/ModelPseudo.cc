/*
 * (C) Copyright 2019-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

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
static oops::interface::ModelMaker<Traits, ModelPseudo> makermodel_("PSEUDO");
// -------------------------------------------------------------------------------------------------
ModelPseudo::ModelPseudo(const Geometry & resol, const Parameters_ & params)
  : tstep_(0),
    io_()
{
  oops::Log::trace() << "ModelPseudo::ModelPseudo starting" << std::endl;
  // Create IO object
  io_.reset(IOFactory::create(resol, *params.ioParametersWrapper.ioParameters.value()));
  // Trace
  // Get timestep from params
  tstep_ = params.tstep.value();
  // Optionally retrieve run stage check
  runstagecheck_ = params.runstagecheck.value();
  // Trace
  oops::Log::trace() << "ModelPseudo::ModelPseudo done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
ModelPseudo::~ModelPseudo() {
  oops::Log::trace() << "ModelPseudo destructed" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelPseudo::initialize(State & xx) const {
  oops::Log::trace() << "ModelPseudo::initialize starting & also done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelPseudo::step(State & xx, const ModelBias &) const {
  xx.validTime() += tstep_;
  if (runstage_) {
    // Read model state at valid time from files
    io_->read(xx);
  } else {
    // Do nothing and print message
    if (oops::mpi::world().rank() == 0) {
      oops::Log::warning() << "Pseudo model has already run through "
                              "once so not re-reading states, just ticking the"
                              " clock." << std::endl;
    }
  }
  oops::Log::trace() << "ModelPseudo::step" << xx.validTime() << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelPseudo::finalize(State & xx) const {
  oops::Log::trace() << "ModelPseudo::finalize starting" << std::endl;
  if (runstagecheck_) {runstage_ = false;}
  oops::Log::trace() << "ModelPseudo::finalize done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelPseudo::print(std::ostream & os) const {
  os << "ModelPseudo::print not implemented";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
