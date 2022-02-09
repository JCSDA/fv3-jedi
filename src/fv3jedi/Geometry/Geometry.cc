/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <mpi.h>

#include <sstream>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/util/Config.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/GeometryIterator/GeometryIterator.interface.h"

// -------------------------------------------------------------------------------------------------

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

Geometry::Geometry(const Parameters_ & params, const eckit::mpi::Comm & comm) :
                   comm_(comm) {
  // Call the initialize phase, done only once.
  static bool initialized = false;
  if (!initialized) {
    fv3jedi_geom_initialize_f90((*params.fmsInit.value()).toConfiguration(), &comm_);
    initialized = true;
    oops::Log::debug() << "FMS MPP initialized on " << comm_.name() << std::endl;
  }

  // Geometry constructor
  int nlev;
  fv3jedi_geom_setup_f90(keyGeom_, params.toConfiguration(), &comm_, nlev);

  // Construct the field sets and add to Geometry
  fieldsMeta_.reset(new FieldsMetadata(params.fieldsMetadataParameters, nlev));
  fv3jedi_geom_addfmd_f90(keyGeom_, fieldsMeta_.get());

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

  // Read the orography
  if (params.orography.value() != boost::none) {
    State orogstate(*this, *params.orography.value());
    orogstate.fillGeomOrography(*this);
  }
}

// -------------------------------------------------------------------------------------------------

Geometry::Geometry(const Geometry & other) : comm_(other.comm_) {
  fieldsMeta_ = std::make_shared<FieldsMetadata>(*other.fieldsMeta_);
  fv3jedi_geom_clone_f90(keyGeom_, other.keyGeom_, fieldsMeta_.get());
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
  int ist, iend, jst, jend, kst, kend, npz;
  fv3jedi_geom_start_end_f90(keyGeom_, ist, iend, jst, jend, kst, kend, npz);
  // 3D iterator starts from 0 for surface variables
  return GeometryIterator(*this, ist, jst, kst);
}

// -------------------------------------------------------------------------------------------------

GeometryIterator Geometry::end() const {
  // return end of the geometry on this mpi tile
  // (returns index out of bounds for the iterator loops to work)

  return GeometryIterator(*this, -1, -1, -1);
}

// -------------------------------------------------------------------------------------------------

std::vector<double> Geometry::verticalCoord(std::string & vcUnits) const {
  // returns vertical coordinate in untis of vcUnits
  // to enable initial comparisons with GSI, verticalCoord is valid for psurf=1e5
  // TODO(TBD) implement interface where vertical coordinate is valid for state[gridIterator].psurf

  int ist, iend, jst, jend, kst, kend, npz;
  fv3jedi_geom_start_end_f90(keyGeom_, ist, iend, jst, jend, kst, kend, npz);
  std::vector<double> vc(npz);
  double psurf = 100000.0;
  if (vcUnits == "logp") {
    fv3jedi_geom_verticalCoord_f90(keyGeom_, vc[0], npz, psurf);
  } else if (vcUnits == "levels") {
    for (int i=0; i < npz; ++i) {vc[i]=i+1;}
  } else {
    std::stringstream errorMsg;
    errorMsg << "Uknown vertical coordinate unit " << vcUnits << std::endl;
    ABORT(errorMsg.str());
  }
  oops::Log::debug() << "fv3 vert coord: " << vc << std::endl;
  return vc;
}

// -------------------------------------------------------------------------------------------------

std::vector<size_t> Geometry::variableSizes(const oops::Variables & vars) const {
  // Array of level heights
  std::vector<size_t> varSizes;
  // Loop through arrays and search metadata map for the levels
  for (size_t it = 0; it < vars.size(); it++) {
    varSizes.push_back(fieldsMeta_->getLevelsFromLongName(vars[it]));
  }
  return varSizes;
}

// -------------------------------------------------------------------------------------------------

void Geometry::print(std::ostream & os) const {
  int cube;
  fv3jedi_geom_print_f90(keyGeom_, cube);
  os << "fv3jedi::Geometry resolution: c" << cube;
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
