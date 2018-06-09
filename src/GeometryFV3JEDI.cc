/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <mpi.h>
#include "src/GeometryFV3JEDI.h"
#include "oops/util/Logger.h"
#include "Fortran.h"
#include "eckit/config/Configuration.h"
#include "UtilitiesFV3JEDI.h"

// -----------------------------------------------------------------------------
namespace fv3jedi {
// -----------------------------------------------------------------------------
GeometryFV3JEDI::GeometryFV3JEDI(const eckit::Configuration & conf) {
  const eckit::Configuration * configc = &conf;
  stageFv3Files(conf);
  fv3jedi_geo_setup_f90(keyGeom_, &configc);
  removeFv3Files();
}
// -----------------------------------------------------------------------------
GeometryFV3JEDI::GeometryFV3JEDI(const GeometryFV3JEDI & other) {
  const int key_geo = other.keyGeom_;
  fv3jedi_geo_clone_f90(key_geo, keyGeom_);
}
// -----------------------------------------------------------------------------
GeometryFV3JEDI::~GeometryFV3JEDI() {
  fv3jedi_geo_delete_f90(keyGeom_);
}
// -----------------------------------------------------------------------------
void GeometryFV3JEDI::print(std::ostream & os) const {
  fv3jedi_geo_info_f90(keyGeom_);
  os << "fv3jedi::GeometryFV3JEDI::print not implemented yet";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
