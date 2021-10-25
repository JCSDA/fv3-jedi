/*
 * (C) Copyright 2019-2021 UCAR
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
#include "fv3jedi/IO/Utils/IOBase.h"
#include "fv3jedi/Model/pseudo/ModelPseudo.h"
#include "fv3jedi/ModelBias/ModelBias.h"
#include "fv3jedi/State/State.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static oops::interface::ModelMaker<Traits, ModelPseudo> makermodel_("PSEUDO");
// -------------------------------------------------------------------------------------------------
ModelPseudo::ModelPseudo(const Geometry & resol, const eckit::Configuration & mconf)
  : tstep_(0), vars_(mconf, "model variables"), io_(IOFactory::create(mconf, resol))
{
  // Trace
  oops::Log::trace() << "ModelPseudo::ModelPseudo" << std::endl;
  // Get timestep from condfig
  tstep_ = util::Duration(mconf.getString("tstep"));
  // Optionally retrieve run stage check
  runstagecheck_ = mconf.getBool("run stage check", false);
  // Trace
  oops::Log::trace() << "ModelPseudo created" << std::endl;
}
// -------------------------------------------------------------------------------------------------
ModelPseudo::~ModelPseudo() {
  oops::Log::trace() << "ModelPseudo destructed" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelPseudo::initialize(State & xx) const {
  oops::Log::trace() << "ModelPseudo::initialize" << std::endl;
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
  if (runstagecheck_) {runstage_ = false;}
  oops::Log::trace() << "ModelPseudo::finalize" << std::endl;
}
// -------------------------------------------------------------------------------------------------
void ModelPseudo::print(std::ostream & os) const {
  os << "ModelPseudo::print not implemented";
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
