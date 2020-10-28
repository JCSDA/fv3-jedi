/*
 * (C) Copyright 2017-2020 UCAR
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

// -------------------------------------------------------------------------------------------------
namespace fv3jedi {
// -------------------------------------------------------------------------------------------------

ErrorCovariance::ErrorCovariance(const Geometry & resol, const oops::Variables &,
                                 const eckit::Configuration & conf, const State &, const State &) {
  oops::Log::trace() << "ErrorCovariance created" << std::endl;
}

// -------------------------------------------------------------------------------------------------

ErrorCovariance::~ErrorCovariance() {
  oops::Log::trace() << "ErrorCovariance destructed" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void ErrorCovariance::multiply(const Increment & dxin, Increment & dxout) const {
  dxout = dxin;  // Identity
}

// -------------------------------------------------------------------------------------------------

void ErrorCovariance::inverseMultiply(const Increment & dxin, Increment & dxout) const {
  dxout = dxin;  // Identity
}

// -------------------------------------------------------------------------------------------------

void ErrorCovariance::randomize(Increment & dx) const {
  dx.random();
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
