/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <ostream>
#include <string>

#include "oops/base/LinearModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/Utilities/Traits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

// Linear model definition.

class Tlm: public oops::LinearModelBase<Traits>,
                private util::ObjectCounter<Tlm> {
 public:
  static const std::string classname() {return "fv3jedi::Tlm";}

  // Constructor/destructor
  Tlm(const Geometry &, const eckit::Configuration &);
  ~Tlm();

  // Set the trajectory
  void setTrajectory(const State &, State &, const ModelBias &) override;

  // Run TLM and its adjoint
  void initializeTL(Increment &) const override;
  void stepTL(Increment &, const ModelBiasIncrement &) const override;
  void finalizeTL(Increment &) const override;

  void initializeAD(Increment &) const override;
  void stepAD(Increment &, ModelBiasIncrement &) const override;
  void finalizeAD(Increment &) const override;

  // Accessor functions
  const util::Duration & timeResolution() const override {return tstep_;}
  const oops::Variables & variables() const override {return linvars_;}

 private:
  void print(std::ostream &) const override;
  typedef std::map< util::DateTime, int >::iterator trajIter;
  typedef std::map< util::DateTime, int >::const_iterator trajICst;

// Data
  F90model keySelf_;
  util::Duration tstep_;
  std::map< util::DateTime, F90traj> trajmap_;
  const oops::Variables linvars_;
};
// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
