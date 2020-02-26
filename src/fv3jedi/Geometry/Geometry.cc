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
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"
#include "atlas/util/Config.h"

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

  // Create ATLAS grid configuration
  const atlas::util::Config atlasConfig;
  const eckit::Configuration * fconf = &atlasConfig;
  fv3jedi_geo_create_atlas_grid_conf_f90(keyGeom_, &fconf);

  // Create ATLAS grid
  atlas::UnstructuredGrid atlasUnstructuredGrid(atlasConfig);

  // Create mesh
  atlas::MeshGenerator atlasMeshGenerator("no_connectivity");
  atlas::Mesh atlasMesh = atlasMeshGenerator.generate(atlasUnstructuredGrid);

  // Create ATLAS function space
  atlasFunctionSpace_.reset(new atlas::functionspace::NodeColumns(atlasMesh,
                            atlas::option::halo(0)));

  // Set ATLAS function space pointer in Fortran
  fv3jedi_geo_set_atlas_functionspace_pointer_f90(keyGeom_, atlasFunctionSpace_->get());

  // Fill ATLAS fieldset
  atlasFieldSet_.reset(new atlas::FieldSet());
  fv3jedi_geo_fill_atlas_fieldset_f90(keyGeom_, atlasFieldSet_->get());
}
// -----------------------------------------------------------------------------
Geometry::Geometry(const Geometry & other) : comm_(other.comm_) {
  const int key_geo = other.keyGeom_;
  fv3jedi_geo_clone_f90(key_geo, keyGeom_);
  atlasFunctionSpace_.reset(new atlas::functionspace::NodeColumns(
                            other.atlasFunctionSpace_->mesh(), atlas::option::halo(0)));
  atlasFieldSet_.reset(new atlas::FieldSet());
  for (int jfield = 0; jfield < other.atlasFieldSet_->size(); ++jfield) {
    atlas::Field atlasField = other.atlasFieldSet_->field(jfield);
    atlasFieldSet_->add(atlasField);
  }
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
