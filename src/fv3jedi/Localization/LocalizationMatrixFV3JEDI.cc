/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Localization/LocalizationMatrixFV3JEDI.h"

#include "eckit/config/Configuration.h"
#include "fv3jedi/Geometry/GeometryFV3JEDI.h"
#include "fv3jedi/Increment/IncrementFV3JEDI.h"
#include "LocalizationMatrixFV3JEDIFortran.h"

// -----------------------------------------------------------------------------
namespace fv3jedi {
// -----------------------------------------------------------------------------
LocalizationMatrixFV3JEDI::LocalizationMatrixFV3JEDI(const GeometryFV3JEDI &
                                                 resol,
                                                 const eckit::Configuration &
                                                 config) {
}
// -----------------------------------------------------------------------------
LocalizationMatrixFV3JEDI::~LocalizationMatrixFV3JEDI() {
}
// -----------------------------------------------------------------------------
void LocalizationMatrixFV3JEDI::multiply(IncrementFV3JEDI & dx) const {
}
// -----------------------------------------------------------------------------
void LocalizationMatrixFV3JEDI::print(std::ostream & os) const {
  os << "LocalizationMatrixFV3JEDI::print not implemented";
}
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
