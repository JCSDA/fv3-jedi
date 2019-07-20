/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_TLM_TLMID_H_
#define FV3JEDI_TLM_TLMID_H_

#include <string>

#include <boost/noncopyable.hpp>

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
/// FV3JEDI linear identity model definition.
/*!
 *  FV3JEDI linear identity model definition and configuration parameters.
 */

class TlmId: public oops::LinearModelBase<Traits>,
                  private util::ObjectCounter<TlmId> {
 public:
  static const std::string classname() {return "fv3jedi::TlmId";}

  TlmId(const Geometry &, const eckit::Configuration &);
  ~TlmId();

/// Model trajectory computation
  void setTrajectory(const State &, State &,
                      const ModelBias &) override;

/// Run TLM and its adjoint
  void initializeTL(Increment &) const override;
  void stepTL(Increment &, const ModelBiasIncrement &)
               const override;
  void finalizeTL(Increment &) const override;

  void initializeAD(Increment &) const override;
  void stepAD(Increment &, ModelBiasIncrement &)
                const override;
  void finalizeAD(Increment &) const override;

/// Other utilities
  const util::Duration & timeResolution() const override {return tstep_;}
  const Geometry & resolution() const {return resol_;}
  const oops::Variables & variables() const override {return linvars_;}

 private:
  void print(std::ostream &) const override;

// Data
  int keyConfig_;
  util::Duration tstep_;
  const Geometry resol_;
  const oops::Variables linvars_;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_TLM_TLMID_H_
