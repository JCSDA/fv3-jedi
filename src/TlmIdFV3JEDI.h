/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef FV3JEDI_MODEL_TLMIDFV3JEDI_H_
#define FV3JEDI_MODEL_TLMIDFV3JEDI_H_

#include <string>

#include <boost/noncopyable.hpp>

#include "oops/interface/LinearModelBase.h"

#include "util/Duration.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"

#include "FV3JEDITraits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {

// -----------------------------------------------------------------------------
/// FV3JEDI linear identity model definition.
/*!
 *  FV3JEDI linear identity model definition and configuration parameters.
 */

class TlmIdFV3JEDI: public oops::LinearModelBase<FV3JEDITraits>,
                  private util::ObjectCounter<TlmIdFV3JEDI> {
 public:
  static const std::string classname() {return "fv3jedi::TlmIdFV3JEDI";}

  TlmIdFV3JEDI(const GeometryFV3JEDI &, const eckit::Configuration &);
  ~TlmIdFV3JEDI();

/// Model trajectory computation
  void setTrajectory(const StateFV3JEDI &, StateFV3JEDI &, const ModelBiasFV3JEDI &) override;

/// Run TLM and its adjoint
  void initializeTL(IncrementFV3JEDI &) const override;
  void stepTL(IncrementFV3JEDI &, const ModelBiasIncrementFV3JEDI &) const override;
  void finalizeTL(IncrementFV3JEDI &) const override;

  void initializeAD(IncrementFV3JEDI &) const override;
  void stepAD(IncrementFV3JEDI &, ModelBiasIncrementFV3JEDI &) const override;
  void finalizeAD(IncrementFV3JEDI &) const override;

/// Other utilities
  const util::Duration & timeResolution() const override {return tstep_;}
  const GeometryFV3JEDI & resolution() const {return resol_;}

 private:
  void print(std::ostream &) const override;

// Data
  int keyConfig_;
  util::Duration tstep_;
  const GeometryFV3JEDI resol_;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_MODEL_TLMIDFV3JEDI_H_
