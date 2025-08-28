/*
 * (C) Copyright 2025-2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

#include "fv3jedi/Geometry/Geometry.h"

namespace fv3jedi {
  class ModelBias;
  class State;

// -----------------------------------------------------------------------------

class ModelBase : public util::Printable {
 public:
  ModelBase() = default;
  virtual ~ModelBase() = default;

  virtual void initialize(State &) const = 0;
  virtual void step(State &, const ModelBias &) const = 0;
  virtual void finalize(State &) const = 0;

  virtual const util::Duration & timeResolution() const = 0;

 private:
  void print(std::ostream &) const override = 0;
};

// -----------------------------------------------------------------------------

class ModelFactory {
 public:
  static ModelBase * create(const Geometry &, const eckit::Configuration &);

  virtual ~ModelFactory() = default;

 protected:
  explicit ModelFactory(const std::string & name);

 private:
  virtual ModelBase * make(const Geometry &, const eckit::Configuration &) = 0;

  static std::map < std::string, ModelFactory * > & getMakers() {
    static std::map < std::string, ModelFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class ModelMaker : public ModelFactory {
 public:
  explicit ModelMaker(const std::string & name) : ModelFactory(name) {}

  ModelBase * make(const Geometry & geom, const eckit::Configuration & config) override {
    oops::Log::trace() << "ModelBase::make starting" << std::endl;
    return new T(geom, config);
  }
};

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
