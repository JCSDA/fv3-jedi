/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <mpi.h>

#include "eckit/config/Configuration.h"
#include "oops/util/Logger.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/GeometryIterator/GeometryIterator.interface.h"
#include "fv3jedi/Run/Run.h"
#include "fv3jedi/Utilities/Utilities.h"

// -----------------------------------------------------------------------------
namespace fv3jedi {
// -----------------------------------------------------------------------------
Geometry::Geometry(const eckit::Configuration & conf,
                   const eckit::mpi::Comm & comm) : comm_(comm) {
  const eckit::Configuration * configc = &conf;

  static bool initialized = false;
  if (!initialized) {
    stageFMSFiles(conf);
    fv3jedi_setup_f(&configc, &comm_);
    removeFv3Files();
    initialized = true;
    oops::Log::debug() << "FMS MPP initialized on " << comm_.name() << std::endl;
  }

  stageFv3Files(conf);
  fv3jedi_geo_setup_f90(keyGeom_, &configc, &comm_);
  removeFv3Files();
}
// -----------------------------------------------------------------------------
Geometry::Geometry(const Geometry & other) : comm_(other.comm_) {
  const int key_geo = other.keyGeom_;
  fv3jedi_geo_clone_f90(key_geo, keyGeom_);
}
// -----------------------------------------------------------------------------
Geometry::~Geometry() {
  fv3jedi_geo_delete_f90(keyGeom_);
}
// -----------------------------------------------------------------------------
GeometryIterator Geometry::begin() const {
  // return start of the geometry on this mpi tile
  int ist, iend, jst, jend, npz;
  fv3jedi_geo_start_end_f90(keyGeom_, ist, iend, jst, jend, npz);
  return GeometryIterator(*this, ist, jst);
}
// -----------------------------------------------------------------------------
GeometryIterator Geometry::end() const {
  // return end of the geometry on this mpi tile
  // (returns index out of bounds for the iterator loops to work)

  return GeometryIterator(*this, -1, -1);
}
// -----------------------------------------------------------------------------
void Geometry::print(std::ostream & os) const {
  fv3jedi_geo_info_f90(keyGeom_);
  os << "fv3jedi::Geometry::print not implemented yet";
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
