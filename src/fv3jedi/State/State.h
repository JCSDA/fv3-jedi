/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/base/WriteParametersBase.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "fv3jedi/IO/Utils/IOBase.h"
#include "fv3jedi/State/State.interface.h"

namespace oops {
  class DateTime;
  class Variables;
}

namespace fv3jedi {
  class Geometry;
  class Increment;

// -------------------------------------------------------------------------------------------------

class AnalyticICParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(AnalyticICParameters, Parameters)
 public:
  // Analytic initial condition parameters
  oops::RequiredParameter<std::string> method{ "method", this};
};

// -------------------------------------------------------------------------------------------------

class StateParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(StateParameters, Parameters)
 public:
  // Analytic initial condition parameters
  oops::OptionalParameter<oops::Variables> stateVariables{ "state variables", this};
  oops::OptionalParameter<AnalyticICParameters> analytic{ "analytic init", this};
  oops::OptionalParameter<util::DateTime> datetime{"datetime", this};
  // Read parameters
  IOParametersWrapper ioParametersWrapper{this};
  oops::OptionalParameter<bool> setdatetime{"set datetime on read", this};
};

// -------------------------------------------------------------------------------------------------

class StateWriteParameters : public oops::WriteParametersBase {
  OOPS_CONCRETE_PARAMETERS(StateWriteParameters, WriteParametersBase)
 public:
  IOParametersWrapper ioParametersWrapper{this};
};

// -------------------------------------------------------------------------------------------------

class State : public util::Printable, private util::ObjectCounter<State> {
 public:
  static const std::string classname() {return "fv3jedi::State";}

  typedef StateParameters Parameters_;
  typedef StateWriteParameters WriteParameters_;
  typedef AnalyticICParameters AnalyticICParameters_;

  typedef std::unique_ptr<IOBase> IOBase_;

// Constructor, destructor and basic operators
  State(const Geometry &, const oops::Variables &, const util::DateTime &);
  State(const Geometry &, const Parameters_ &);
  State(const Geometry &, const State &);
  State(const State &);
  virtual ~State();

  State & operator=(const State &);
  void zero();
  void accumul(const double &, const State &);

// Interpolate state
  void changeResolution(const State & xx);

// Interactions with Increment
  State & operator+=(const Increment &);

// IO and diagnostics
  void analytic_init(const AnalyticICParameters_ &, const Geometry &);
  void fillGeomOrography(Geometry &) const;
  void read(const Parameters_ &);
  void write(const WriteParameters_ &) const;
  double norm() const;

/// Serialize and deserialize
  size_t serialSize() const;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, size_t &);

// Add or remove fields
  void updateFields(const oops::Variables &);

// Utilities
  std::shared_ptr<const Geometry> geometry() const {return geom_;}
  const oops::Variables & variables() const {return vars_;}
  const oops::Variables & variablesLongName() const {return varsLongName_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}
  void updateTime(const util::Duration & dt) {time_ += dt;}

  int & toFortran() {return keyState_;}
  const int & toFortran() const {return keyState_;}

// Private methods and variables
 private:
  void print(std::ostream &) const;
  F90state keyState_;
  std::shared_ptr<const Geometry> geom_;
  oops::Variables vars_;
  oops::Variables varsLongName_;
  util::DateTime time_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
