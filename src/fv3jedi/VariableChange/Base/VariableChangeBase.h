/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/VariableChangeParametersBase.h"
#include "oops/base/Variables.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParametersOrConfiguration.h"
#include "oops/util/parameters/PolymorphicParameter.h"
#include "oops/util/Printable.h"
#include "vader/VaderParameters.h"

namespace fv3jedi {
  class Geometry;
  class State;

// -------------------------------------------------------------------------------------------------

class VariableChangeParametersBase : public oops::VariableChangeParametersBase {
  OOPS_ABSTRACT_PARAMETERS(VariableChangeParametersBase, oops::VariableChangeParametersBase)

 public:
  oops::OptionalParameter<std::string> name{"variable change name", this};
  oops::Parameter<std::map<std::string, std::vector<std::string>>> vaderCustomCookbook{
    "vader custom cookbook",
    // Default VADER cookbook for fv3-jedi
    {
      // pt: from t and pkz
      {"potential_temperature",        {"AirPotentialTemperature_B"}},
      // P: from delp, from ps (and ak/bk)
      {"air_pressure_levels",          {"AirPressureAtInterface_B", "AirPressureAtInterface_A"}},
      // p: from pe
      {"air_pressure",                 {"AirPressure_A"}},
      // ln(p) from pe
      {"ln_air_pressure_at_interface", {"LnAirPressureAtInterface_A"}},
      // p^kappa from pe and ln(p)
      {"air_pressure_to_kappa",        {"AirPressureToKappa_A"}},
      // delp: from p
      {"air_pressure_thickness",       {"AirPressureThickness_A"}},
      // pt: from t and ps
      {"potential_temperature",        {"AirPotentialTemperature_A"}},
      // ps: from delp
      {"surface_pressure",             {"SurfaceAirPressure_A"}},
      // tv: from t and q
      {"virtual_temperature",          {"AirVirtualTemperature_A"}}
    },
    this};
  oops::Parameter<vader::VaderParameters> vader{"vader", {}, this};
};

// -------------------------------------------------------------------------------------------------

class GenericVariableChangeParameters : public VariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericVariableChangeParameters, VariableChangeParametersBase)
 public:
  oops::ConfigurationParameter config{this};
};

// -------------------------------------------------------------------------------------------------

class VariableChangeBase : public util::Printable, private boost::noncopyable {
 public:
  VariableChangeBase() {}
  virtual ~VariableChangeBase() {}

  virtual void changeVar(const State &, State &) const = 0;
  virtual void changeVarInverse(const State &, State &) const = 0;

 private:
  virtual void print(std::ostream &) const = 0;
};

// -------------------------------------------------------------------------------------------------

class VariableChangeFactory;

// -------------------------------------------------------------------------------------------------

class VariableChangeParametersWrapper : public oops::VariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(VariableChangeParametersWrapper, oops::VariableChangeParametersBase)
 public:
  oops::PolymorphicParameter<fv3jedi::VariableChangeParametersBase, VariableChangeFactory>
    variableChangeParameters{"variable change name", "default", this};
  // During the transition to Vader it is useful to run with either just vader or just the fv3-jedi
  // variable transforms to avoid not knowing which part of the code is doing the transforms
  oops::Parameter<bool> run_vader{"run vader", true, this};
  oops::Parameter<bool> run_fv3jedi{"run fv3jedi", true, this};
};

// -------------------------------------------------------------------------------------------------

class VariableChangeFactory {
 public:
  static VariableChangeBase * create(const Geometry & geom,
                                      const VariableChangeParametersBase & params);

  static std::unique_ptr<VariableChangeParametersBase> createParameters(const std::string &name);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~VariableChangeFactory() = default;

 protected:
  explicit VariableChangeFactory(const std::string &name);

 private:
  virtual VariableChangeBase * make(const Geometry &, const VariableChangeParametersBase &) = 0;

  virtual std::unique_ptr<VariableChangeParametersBase> makeParameters() const = 0;

  static std::map < std::string, VariableChangeFactory * > & getMakers() {
    static std::map < std::string, VariableChangeFactory * > makers_;
    return makers_;
  }
};

// -------------------------------------------------------------------------------------------------

template<class T>
class VariableChangeMaker : public VariableChangeFactory {
  typedef oops::TParameters_IfAvailableElseFallbackType_t<T, GenericVariableChangeParameters>
    Parameters_;

  VariableChangeBase * make(const Geometry & geom,
                             const VariableChangeParametersBase & params) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return new T(geom, oops::parametersOrConfiguration<oops::HasParameters_<T>::value>(
                   stronglyTypedParams));
  }

  std::unique_ptr<VariableChangeParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit VariableChangeMaker(const std::string & name) : VariableChangeFactory(name) {}
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
