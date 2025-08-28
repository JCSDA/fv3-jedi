/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/LinearVariableChange/Base/LinearVariableChangeBase.h"
#include "fv3jedi/Utilities/interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {
  class ModelBias;
  class ModelBiasIncrement;

// -------------------------------------------------------------------------------------------------

// Linear model definition.

class Tlm: public util::Printable,
           private util::ObjectCounter<Tlm> {
 public:
  static const std::string classname() {return "fv3jedi::Tlm";}
  static std::vector<std::string> names() {return {"FV3JEDITLM"};}

  // Constructor/destructor
  Tlm(const Geometry &, const eckit::Configuration &);
  ~Tlm();

  // Set the trajectory
  void setTrajectory(const State &, State &, const ModelBias &);

  // Run TLM and its adjoint
  void initializeTL(Increment &) const;
  void stepTL(Increment &, const ModelBiasIncrement &) const;
  void finalizeTL(Increment &) const;

  void initializeAD(Increment &) const;
  void stepAD(Increment &, ModelBiasIncrement &) const;
  void finalizeAD(Increment &) const;

  // Accessor functions
  const util::Duration & timeResolution() const {return tstep_;}
  const util::Duration & stepTrajectory() const {return tstep_;}

 private:
  void print(std::ostream &) const override;
  typedef std::map< util::DateTime, int >::iterator trajIter;
  typedef std::map< util::DateTime, int >::const_iterator trajICst;

// Data
  const Geometry & geom_;
  F90model keySelf_;
  util::Duration tstep_;
  std::map<util::DateTime, F90traj> trajmap_;
};
// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
