/*
 * (C) Copyright 2025-2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <memory>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Model/ModelWrapper.h"

namespace fv3jedi {
  class ModelBias;
  class State;

// -----------------------------------------------------------------------------

ModelWrapper::ModelWrapper(const Geometry & resol, const eckit::Configuration & conf) : model_() {
  oops::Log::trace() << "ModelWrapper::Model starting" << std::endl;
  model_.reset(ModelFactory::create(resol, conf));
  oops::Log::trace() << "ModelWrapper::Model done" << std::endl;
}

// -----------------------------------------------------------------------------

void ModelWrapper::initialize(State & xx) const {
  oops::Log::trace() << "ModelWrapper::initialize starting" << std::endl;
  model_->initialize(xx);
  oops::Log::trace() << "ModelWrapper::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

void ModelWrapper::step(State & xx, const ModelBias & maux) const {
  oops::Log::trace() << "ModelWrapper::step starting" << std::endl;
  model_->step(xx, maux);
  oops::Log::trace() << "ModelWrapper::step done" << std::endl;
}

// -----------------------------------------------------------------------------

void ModelWrapper::finalize(State & xx) const {
  oops::Log::trace() << "ModelWrapper::finalize starting" << std::endl;
  model_->finalize(xx);
  oops::Log::trace() << "ModelWrapper::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

void ModelWrapper::print(std::ostream & os) const {
  oops::Log::trace() << "ModelWrapper::print starting" << std::endl;
  os << *model_;
  oops::Log::trace() << "ModelWrapper::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
