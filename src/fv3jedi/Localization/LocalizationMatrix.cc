/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Localization/LocalizationMatrix.h"

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "LocalizationMatrixFortran.h"

// -----------------------------------------------------------------------------
namespace fv3jedi {
// -----------------------------------------------------------------------------
LocalizationMatrix::LocalizationMatrix(const Geometry &
                                                 resol,
                                                 const eckit::Configuration &
                                                 config) {
}
// -----------------------------------------------------------------------------
LocalizationMatrix::~LocalizationMatrix() {
}
// -----------------------------------------------------------------------------
void LocalizationMatrix::multiply(Increment & dx) const {
}
// -----------------------------------------------------------------------------
void LocalizationMatrix::print(std::ostream & os) const {
  os << "LocalizationMatrix::print not implemented";
}
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
