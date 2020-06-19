/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <mpi.h>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/util/Config.h"

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/GeometryIterator/GeometryIterator.interface.h"
#include "fv3jedi/Run/Run.h"
#include "fv3jedi/Utilities/Utilities.h"

// -------------------------------------------------------------------------------------------------

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

Geometry::Geometry(const eckit::Configuration & conf,
                   const eckit::mpi::Comm & comm) : comm_(comm), fieldsMeta_(conf) {
  const eckit::Configuration * configc = &conf;

  static bool initialized = false;
  if (!initialized) {
    stageFMSFiles(conf, comm);
    fv3jedi_setup_f(&configc, &comm_);
    removeFv3Files(comm);
    initialized = true;
    oops::Log::debug() << "FMS MPP initialized on " << comm_.name() << std::endl;
  }

  stageFv3Files(conf, comm);
  if ( !conf.has("nml_file") ) {
    generateGeomFv3Conf(conf, comm);
  }
  fv3jedi_geom_setup_f90(keyGeom_, &configc, &comm_, &fieldsMeta_);
  removeFv3Files(comm);

  // Set ATLAS lon/lat field
  atlasFieldSet_.reset(new atlas::FieldSet());
  fv3jedi_geom_set_atlas_lonlat_f90(keyGeom_, atlasFieldSet_->get());
  atlas::Field atlasField = atlasFieldSet_->field("lonlat");

  // Create ATLAS function space
  atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(atlasField));

  // Set ATLAS function space pointer in Fortran
  fv3jedi_geom_set_atlas_functionspace_pointer_f90(keyGeom_, atlasFunctionSpace_->get());

  // Fill ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());
  fv3jedi_geom_fill_atlas_fieldset_f90(keyGeom_, atlasFieldSet_->get());
}

// -------------------------------------------------------------------------------------------------

Geometry::Geometry(const Geometry & other) : comm_(other.comm_), fieldsMeta_(other.fieldsMeta_) {
  fv3jedi_geom_clone_f90(keyGeom_, other.keyGeom_, &fieldsMeta_);
  atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(
                            other.atlasFunctionSpace_->lonlat()));
  fv3jedi_geom_set_atlas_functionspace_pointer_f90(keyGeom_, atlasFunctionSpace_->get());
  atlasFieldSet_.reset(new atlas::FieldSet());
  for (int jfield = 0; jfield < other.atlasFieldSet_->size(); ++jfield) {
    atlas::Field atlasField = other.atlasFieldSet_->field(jfield);
    atlasFieldSet_->add(atlasField);
  }
}

// -------------------------------------------------------------------------------------------------

Geometry::~Geometry() {
  fv3jedi_geom_delete_f90(keyGeom_);
}

// -------------------------------------------------------------------------------------------------

GeometryIterator Geometry::begin() const {
  // return start of the geometry on this mpi tile
  int ist, iend, jst, jend, npz;
  fv3jedi_geom_start_end_f90(keyGeom_, ist, iend, jst, jend, npz);
  return GeometryIterator(*this, ist, jst);
}

// -------------------------------------------------------------------------------------------------

GeometryIterator Geometry::end() const {
  // return end of the geometry on this mpi tile
  // (returns index out of bounds for the iterator loops to work)

  return GeometryIterator(*this, -1, -1);
}

// -------------------------------------------------------------------------------------------------

void Geometry::print(std::ostream & os) const {
  int cube;
  fv3jedi_geom_print_f90(keyGeom_, cube);
  os << "fv3jedi::Geometry resolution: c" << cube;
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
