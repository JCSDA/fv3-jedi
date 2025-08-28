/*
 * (C) Copyright 2025-2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "oops/util/Printable.h"

#include "fv3jedi/Model/ModelBase.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class Geometry;
  class ModelBias;
  class State;

// -----------------------------------------------------------------------------

class ModelWrapper : public util::Printable {
 public:
  static std::vector<std::string> names() {return {"FV3LM", "GEOS", "PSEUDO", "UFS"};}

  ModelWrapper(const Geometry &, const eckit::Configuration &);
  ~ModelWrapper() = default;

  void initialize(State &) const;
  void step(State &, const ModelBias &) const;
  void finalize(State &) const;

  const util::Duration & timeResolution() const {return model_->timeResolution();}

 private:
  void print(std::ostream &) const override;

  std::unique_ptr<ModelBase> model_;
};

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
