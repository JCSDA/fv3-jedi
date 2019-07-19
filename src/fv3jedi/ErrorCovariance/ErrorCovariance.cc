/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <cmath>

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "fv3jedi/ErrorCovariance/ErrorCovariance.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"

// -----------------------------------------------------------------------------
namespace fv3jedi {
// -----------------------------------------------------------------------------

ErrorCovariance::ErrorCovariance(const Geometry & resol,
                                     const oops::Variables &,
                                     const eckit::Configuration & conf,
                                     const State &,
                                     const State &) {
  time_ = util::DateTime(conf.getString("date"));
  const eckit::Configuration * configc = &conf;
  fv3jedi_b_setup_f90(keyFtnConfig_, &configc, resol.toFortran());
  oops::Log::trace() << "ErrorCovariance created" << std::endl;
}

// -----------------------------------------------------------------------------

ErrorCovariance::~ErrorCovariance() {
  fv3jedi_b_delete_f90(keyFtnConfig_);
  oops::Log::trace() << "ErrorCovariance destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ErrorCovariance::linearize(const State &,
                                       const Geometry & resol) {
  geom_.reset(new Geometry(resol));
}

// -----------------------------------------------------------------------------

void ErrorCovariance::multiply(const Increment & dxin,
                                      Increment & dxout) const {
    dxout = dxin;
//  fv3jedi_b_mult_f90(keyFtnConfig_, dxin.toFortran(),
//                            dxout.toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovariance::inverseMultiply(const Increment & dxin,
                                           Increment & dxout) const {
    dxout = dxin;
//  fv3jedi_b_invmult_f90(keyFtnConfig_, dxin.toFortran(),
//                               dxout.toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovariance::randomize(Increment & dx) const {
  fv3jedi_b_randomize_f90(keyFtnConfig_, dx.toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovariance::print(std::ostream & os) const {
  os << "ErrorCovariance::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
