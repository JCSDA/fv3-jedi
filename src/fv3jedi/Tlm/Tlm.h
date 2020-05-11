/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_TLM_TLM_H_
#define FV3JEDI_TLM_TLM_H_

#include <map>
#include <ostream>
#include <string>

#include "oops/base/LinearModelBase.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/Tlm/Tlm.interface.h"
#include "fv3jedi/Utilities/Traits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {

// -----------------------------------------------------------------------------
/// FV3JEDI linear model definition.
/*!
 *  FV3JEDI linear model definition and configuration parameters.
 */

class Tlm: public oops::LinearModelBase<Traits>,
                private util::ObjectCounter<Tlm> {
 public:
  static const std::string classname() {return "fv3jedi::Tlm";}

  Tlm(const Geometry &, const eckit::Configuration &);
  ~Tlm();

/// Model trajectory computation
  void setTrajectory(const State &, State &,
                     const ModelBias &) override;

/// Run TLM and its adjoint
  void initializeTL(Increment &) const override;
  void stepTL(Increment &, const ModelBiasIncrement &)
               const override;
  void finalizeTL(Increment &) const override;

  void initializeAD(Increment &) const override;
  void stepAD(Increment &, ModelBiasIncrement &) const override;
  void finalizeAD(Increment &) const override;

/// Other utilities
  const util::Duration & timeResolution() const override {return tstep_;}
  const Geometry & resolution() const {return resol_;}
  const oops::Variables & variables() const override {return linvars_;}

 private:
  void print(std::ostream &) const override;
  typedef std::map< util::DateTime, int >::iterator trajIter;
  typedef std::map< util::DateTime, int >::const_iterator trajICst;

// Data
  F90model keyConfig_;
  util::Duration tstep_;
  const Geometry resol_;
  std::map< util::DateTime, F90traj> traj_;
  const ModelTraj lrmodel_;
  const oops::Variables linvars_;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_TLM_TLM_H_
