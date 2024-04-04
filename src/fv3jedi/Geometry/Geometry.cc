/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"

#include "oops/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Geometry/GeometryParameters.h"
#include "fv3jedi/Geometry/TimeInvariantFieldsHelpers.h"
#include "fv3jedi/GeometryIterator/GeometryIterator.interface.h"

// -------------------------------------------------------------------------------------------------

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

Geometry::Geometry(const eckit::Configuration & config, const eckit::mpi::Comm & comm) :
                   comm_(comm), ak_(), bk_() {
  GeometryParameters params;
  params.deserialize(config);
  // Call the initialize phase, done only once.
  static bool initialized = false;
  if (!initialized) {
    fv3jedi_geom_initialize_f90((*params.fmsInit.value()).toConfiguration(), &comm_);
    initialized = true;
  }

  // Geometry constructor
  fv3jedi_geom_setup_f90(keyGeom_, params.toConfiguration(), &comm_, nLevels_);

  // Construct the field sets and add to Geometry
  fieldsMeta_.reset(new FieldsMetadata(params.fieldsMetadataParameters, nLevels_));
  fv3jedi_geom_addfmd_f90(keyGeom_, fieldsMeta_.get());

  {
    // grab local coordinates + global indices from fortran (owned + 1-deep halo points)
    int num_nodes;
    int num_tri_elements;
    int num_quad_elements;
    fv3jedi_geom_get_num_nodes_and_elements_f90(keyGeom_, num_nodes,
                                                num_tri_elements, num_quad_elements);

    std::vector<double> lons(num_nodes);
    std::vector<double> lats(num_nodes);
    std::vector<int> ghosts(num_nodes);
    std::vector<int> global_indices(num_nodes);
    std::vector<int> remote_indices(num_nodes);
    std::vector<int> partitions(num_nodes);

    const int num_tri_boundary_nodes = 3*num_tri_elements;
    const int num_quad_boundary_nodes = 4*num_quad_elements;
    std::vector<int> raw_tri_boundary_nodes(num_tri_boundary_nodes);
    std::vector<int> raw_quad_boundary_nodes(num_quad_boundary_nodes);
    fv3jedi_geom_get_coords_and_connectivities_f90(keyGeom_,
        num_nodes, lons.data(), lats.data(),
        ghosts.data(), global_indices.data(), remote_indices.data(), partitions.data(),
        num_tri_boundary_nodes, raw_tri_boundary_nodes.data(),
        num_quad_boundary_nodes, raw_quad_boundary_nodes.data());

    const int num_elements = num_tri_elements + num_quad_elements;
    std::vector<int> num_elements_per_rank(comm_.size());
    comm_.allGather(num_elements, num_elements_per_rank.begin(), num_elements_per_rank.end());
    int global_element_index = 0;
    for (size_t i = 0; i < comm_.rank(); ++i) {
      global_element_index += num_elements_per_rank[i];
    }

    using atlas::gidx_t;
    using atlas::idx_t;

    std::vector<std::array<gidx_t, 3>> tri_boundary_nodes(num_tri_elements);
    std::vector<gidx_t> tri_global_indices(num_tri_elements);
    for (size_t tri = 0; tri < num_tri_elements; ++tri) {
      for (size_t i = 0; i < 3; ++i) {
        tri_boundary_nodes[tri][i] = raw_tri_boundary_nodes[3*tri + i];
      }
      tri_global_indices[tri] = global_element_index;
      ++global_element_index;
    }
    std::vector<std::array<gidx_t, 4>> quad_boundary_nodes(num_quad_elements);
    std::vector<gidx_t> quad_global_indices(num_quad_elements);
    for (size_t quad = 0; quad < num_quad_elements; ++quad) {
      for (size_t i = 0; i < 4; ++i) {
        quad_boundary_nodes[quad][i] = raw_quad_boundary_nodes[4*quad + i];
      }
      quad_global_indices[quad] = global_element_index;
      ++global_element_index;
    }

    std::vector<atlas::gidx_t> atlas_global_indices(num_nodes);
    std::transform(global_indices.begin(), global_indices.end(), atlas_global_indices.begin(),
                   [](const int index) {return atlas::gidx_t{index};});

    const atlas::idx_t remote_index_base = 1;  // 1-based indexing from Fortran
    std::vector<atlas::idx_t> atlas_remote_indices(num_nodes);
    std::transform(remote_indices.begin(), remote_indices.end(), atlas_remote_indices.begin(),
                   [](const int index) {return atlas::idx_t{index};});

    eckit::LocalConfiguration config{};
    config.set("mpi_comm", comm_.name());

    // establish connectivity
    const atlas::mesh::MeshBuilder mesh_builder{};
    atlas::Mesh mesh = mesh_builder(
        lons,
        lats,
        ghosts,
        atlas_global_indices,
        atlas_remote_indices,
        remote_index_base,
        partitions,
        tri_boundary_nodes,
        tri_global_indices,
        quad_boundary_nodes,
        quad_global_indices,
        config);

    atlas::mesh::actions::build_halo(mesh, 1);
    functionSpace_ = atlas::functionspace::NodeColumns(mesh, config);
  }

  // Set lon/lat field, include halo so we can set up function spaces with/without halo
  atlas::FieldSet fs;
  const bool include_halo = true;
  fv3jedi_geom_set_lonlat_f90(keyGeom_, fs.get(), include_halo);

  // Create function space without halo, for constructing the bump interpolator from fv3jedi
  const atlas::Field lonlatFieldForBump = fs.field("lonlat");
  functionSpaceForBump_ = atlas::functionspace::PointCloud(lonlatFieldForBump);

  // Set function space pointers in Fortran
  fv3jedi_geom_set_functionspace_pointer_f90(keyGeom_, functionSpace_.get(),
                                             functionSpaceForBump_.get());

  // Fill geometry fields. This contains both SABER-related fields and any fields requested to be
  // read from state files in the yamls.
  fields_ = atlas::FieldSet();
  // Add SABER fields
  if (params.timeInvariantFields.value() != boost::none) {
    const auto & timeInvFieldsParams = params.timeInvariantFields.value().value();
    State timeInvState(*this, timeInvFieldsParams.stateFields.value().toConfiguration());
    // Add fields read directly from file
    atlas::FieldSet timeInvFieldSet{};
    timeInvState.toFieldSet(timeInvFieldSet);
    // Populate the atlas halos of the background time-invariant fields. This allows these fields
    // to be used in setting up halo values of derived fields (e.g. the halo of the JEDI sea mask
    // depends on the halo of the UFS slmsk + sheleg).
    timeInvFieldSet.haloExchange();
    for (const auto & f : timeInvFieldSet) {
      fields_.add(f);
    }
    // Compute and add derived time-invariant fields
    if (timeInvFieldsParams.derivedFields.value() != boost::none) {
      const oops::Variables derivedFields = timeInvFieldsParams.derivedFields.value().value();
      insertDerivedTimeInvariantFields(fields_, derivedFields);
    }
  }
  fv3jedi_geom_set_and_fill_geometry_fields_f90(keyGeom_, fields_.get());

  // Copy some Fortran data to C++
  ak_.resize(nLevels_+1);
  bk_.resize(nLevels_+1);
  fv3jedi_geom_get_data_f90(keyGeom_, nLevels_, ak_.data(), bk_.data(), pTop_);
}

// -------------------------------------------------------------------------------------------------

Geometry::Geometry(const Geometry & other) : comm_(other.comm_), ak_(other.ak_), bk_(other.bk_),
nLevels_(other.nLevels_), pTop_(other.pTop_) {
  fieldsMeta_ = std::make_shared<FieldsMetadata>(*other.fieldsMeta_);
  fv3jedi_geom_clone_f90(keyGeom_, other.keyGeom_, fieldsMeta_.get());
  functionSpace_ = atlas::functionspace::NodeColumns(other.functionSpace_);
  functionSpaceForBump_ = atlas::functionspace::PointCloud(other.functionSpaceForBump_.lonlat());
  fv3jedi_geom_set_functionspace_pointer_f90(keyGeom_, functionSpace_.get(),
                                             functionSpaceForBump_.get());
  fields_ = atlas::FieldSet();
  for (auto & field : other.fields_) {
    fields_->add(field);
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
  return vc;
}

// -------------------------------------------------------------------------------------------------

std::vector<size_t> Geometry::variableSizes(const oops::Variables & vars) const {
  // Array of level heights
  std::vector<size_t> varSizes;
  // Loop through arrays and search metadata map for the levels
  for (size_t it = 0; it < vars.size(); it++) {
    varSizes.push_back(fieldsMeta_->getLevels(vars[it]));
  }
  return varSizes;
}

// -------------------------------------------------------------------------------------------------

void Geometry::latlon(std::vector<double> & lats, std::vector<double> & lons,
                      const bool halo) const {
  const atlas::FunctionSpace * fspace;
  if (halo) {
    fspace = &functionSpace_;
  } else {
    fspace = &functionSpaceForBump_;
  }
  const auto lonlat = atlas::array::make_view<double, 2>(fspace->lonlat());
  const size_t npts = fspace->size();
  lats.resize(npts);
  lons.resize(npts);
  for (size_t jj = 0; jj < npts; ++jj) {
    lats[jj] = lonlat(jj, 1);
    lons[jj] = lonlat(jj, 0);
    if (lons[jj] < 0.0) lons[jj] += 360.0;
  }
}

// -------------------------------------------------------------------------------------------------

void Geometry::print(std::ostream & os) const {
  int cube;
  fv3jedi_geom_print_f90(keyGeom_, cube);
  os << "fv3jedi::Geometry resolution: c" << cube;
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
