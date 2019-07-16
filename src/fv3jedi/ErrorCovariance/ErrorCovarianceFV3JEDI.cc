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

#include "fv3jedi/ErrorCovariance/ErrorCovarianceFV3JEDI.h"
#include "fv3jedi/Geometry/GeometryFV3JEDI.h"
#include "fv3jedi/Increment/IncrementFV3JEDI.h"
#include "fv3jedi/State/StateFV3JEDI.h"
#include "ErrorCovarianceFV3JEDIFortran.h"

// -----------------------------------------------------------------------------
namespace fv3jedi {
// -----------------------------------------------------------------------------

ErrorCovarianceFV3JEDI::ErrorCovarianceFV3JEDI(const GeometryFV3JEDI & resol,
                                     const oops::Variables &,
                                     const eckit::Configuration & conf,
                                     const StateFV3JEDI &,
                                     const StateFV3JEDI &) {
  time_ = util::DateTime(conf.getString("date"));
  const eckit::Configuration * configc = &conf;
  fv3jedi_b_setup_f90(keyFtnConfig_, &configc, resol.toFortran());
  oops::Log::trace() << "ErrorCovarianceFV3JEDI created" << std::endl;
}

// -----------------------------------------------------------------------------

ErrorCovarianceFV3JEDI::~ErrorCovarianceFV3JEDI() {
  fv3jedi_b_delete_f90(keyFtnConfig_);
  oops::Log::trace() << "ErrorCovarianceFV3JEDI destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ErrorCovarianceFV3JEDI::linearize(const StateFV3JEDI &,
                                       const GeometryFV3JEDI & resol) {
  geom_.reset(new GeometryFV3JEDI(resol));
}

// -----------------------------------------------------------------------------

void ErrorCovarianceFV3JEDI::multiply(const IncrementFV3JEDI & dxin,
                                      IncrementFV3JEDI & dxout) const {
    dxout = dxin;
//  fv3jedi_b_mult_f90(keyFtnConfig_, dxin.toFortran(),
//                            dxout.toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovarianceFV3JEDI::inverseMultiply(const IncrementFV3JEDI & dxin,
                                           IncrementFV3JEDI & dxout) const {
    dxout = dxin;
//  fv3jedi_b_invmult_f90(keyFtnConfig_, dxin.toFortran(),
//                               dxout.toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovarianceFV3JEDI::randomize(IncrementFV3JEDI & dx) const {
  fv3jedi_b_randomize_f90(keyFtnConfig_, dx.toFortran());
}

// -----------------------------------------------------------------------------

void ErrorCovarianceFV3JEDI::print(std::ostream & os) const {
  os << "ErrorCovarianceFV3JEDI::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
