/*
 * (C) Copyright 2017-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

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

  typedef std::unique_ptr<IOBase> IOBase_;

// Constructor, destructor and basic operators
  State(const Geometry &, const oops::Variables &, const util::DateTime &);
  State(const Geometry &, const eckit::Configuration &);
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
  void analytic_init(const eckit::Configuration &, const Geometry &);
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

// Serialize and deserialize
  size_t serialSize() const;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, size_t &);

// Add or remove fields
  void updateFields(const oops::Variables &);

// Utilities
  const Geometry & geometry() const {return geom_;}
  const oops::Variables & variables() const {return varsJedi_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}
  void updateTime(const util::Duration & dt) {time_ += dt;}

// Accessors to the ATLAS fieldset
  void toFieldSet(atlas::FieldSet &) const;
  void fromFieldSet(const atlas::FieldSet &);

  int & toFortran() {return keyState_;}
  const int & toFortran() const {return keyState_;}

  // Const w.r.t. JEDI, but does update internal fortran state (i.e., the interface-specific fields)
  // to synchronize it with the JEDI-presented fields.
  void synchronizeInterfaceFields() const;
  void setInterfaceFieldsOutOfDate(bool) const;
  const oops::Variables & variablesIncludingInterfaceFields() const {return vars_;}

// Private methods and variables
 private:
  void print(std::ostream &) const;
  F90state keyState_;
  const Geometry & geom_;
  oops::Variables vars_;
  oops::Variables varsJedi_;  // subset of vars_; excluding interface-specific variables
  util::DateTime time_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
