/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <netcdf.h>

#include <algorithm>
#include <chrono>  // NOLINT(build/c++11)
#include <cmath>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "oops/base/GeometryData.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/stringFunctions.h"
#include "oops/util/Timer.h"

#include "fv3jedi/FieldMetadata/FieldsMetadata.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/IO/StructuredGrid/IOStructuredGrid.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Utilities/fv3jedi_vertical_remap.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static IOMaker<IOStructuredGrid> makerIOStructuredGrid_("structured grid");
static IOMaker<IOStructuredGrid> makerIOAuxGrid_("auxgrid");
// -------------------------------------------------------------------------------------------------
static inline void nc_rc(const int return_code, const std::string & operation) {
  if (return_code != NC_NOERR) {
    ABORT("IOStructuredGrid netCDF operation \'" + operation + "\' failed with error: "
          + nc_strerror(return_code));
  }
}

static inline bool validMode(const std::string & mode) {
  return mode == "read" || mode == "write";
}

static inline std::string validModeMessage(const std::string & mode) {
  return "IOStructuredGrid: invalid mode '" + mode
       + "'. Expected one of: 'read' or 'write'.";
}

struct RegionalSubset {
  size_t lat_start;
  size_t lat_count;
  size_t lon_start;
  size_t lon_count;
};

static RegionalSubset computeRegionalSubset(const std::vector<double> & file_lats,
                                            const std::vector<double> & file_lons,
                                            const size_t nLat, const size_t nLon,
                                            const bool fileNorthToSouth,
                                            const double lat_lo_user,
                                            const double lat_hi_user,
                                            const double lon_lo_user,
                                            const double lon_hi_user,
                                            const size_t halo) {
  auto bound_index = [](const std::vector<double> & a, double v, bool left) -> size_t {
    const size_t n = a.size();
    if (n == 0) return 0;
    if (left) {
      if (v <= a.front()) return 0;
      if (v >= a.back()) return n - 1;
      auto it = std::lower_bound(a.begin(), a.end(), v);
      return static_cast<size_t>(std::distance(a.begin(), it));
    }
    if (v <= a.front()) return 0;
    if (v >= a.back()) return n - 1;
    auto it = std::upper_bound(a.begin(), a.end(), v);
    if (it == a.begin()) return 0;
    --it;
    return static_cast<size_t>(std::distance(a.begin(), it));
  };

  double search_lon_min = lon_lo_user;
  double search_lon_max = lon_hi_user;
  if (!file_lons.empty() && file_lons.back() > 180.0) {
    oops::Log::info() << "User bounds: lon=[" << lon_lo_user << "," << lon_hi_user
                      << "] lat=[" << lat_lo_user << "," << lat_hi_user << "]" << std::endl;
    if (search_lon_min < 0.0) search_lon_min += 360.0;
    if (search_lon_max < 0.0) search_lon_max += 360.0;
    oops::Log::info() << "Converted to file lon=[0..360): lon=[" << search_lon_min
                      << "," << search_lon_max << "]" << std::endl;
    if (search_lon_max < search_lon_min) {
      ABORT("Dateline crossing subset is not supported by current subset logic");
    }
  }

  std::vector<double> lats_asc = file_lats;
  if (fileNorthToSouth) std::reverse(lats_asc.begin(), lats_asc.end());

  size_t idx_lat_lo = bound_index(lats_asc, lat_lo_user, true);
  size_t idx_lat_hi = bound_index(lats_asc, lat_hi_user, false);
  if (idx_lat_lo > idx_lat_hi) std::swap(idx_lat_lo, idx_lat_hi);

  size_t idx_lon_lo = bound_index(file_lons, search_lon_min, true);
  size_t idx_lon_hi = bound_index(file_lons, search_lon_max, false);
  if (idx_lon_lo > idx_lon_hi) std::swap(idx_lon_lo, idx_lon_hi);

  idx_lat_lo = (idx_lat_lo > halo) ? (idx_lat_lo - halo) : 0;
  idx_lat_hi = std::min(idx_lat_hi + halo, nLat - 1);
  idx_lon_lo = (idx_lon_lo > halo) ? (idx_lon_lo - halo) : 0;
  idx_lon_hi = std::min(idx_lon_hi + halo, nLon - 1);

  size_t lat_start_file = idx_lat_lo;
  size_t lat_end_file = idx_lat_hi;
  if (fileNorthToSouth) {
    lat_start_file = (nLat - 1) - idx_lat_hi;
    lat_end_file = (nLat - 1) - idx_lat_lo;
  }

  RegionalSubset subset;
  subset.lat_start = lat_start_file;
  subset.lat_count = lat_end_file - lat_start_file + 1;
  subset.lon_start = idx_lon_lo;
  subset.lon_count = idx_lon_hi - idx_lon_lo + 1;

  const size_t lat_end_print = subset.lat_start + subset.lat_count - 1;
  const size_t lon_end_print = subset.lon_start + subset.lon_count - 1;
  oops::Log::info()
    << "Regional subset indices: "
    << "  Latitude:  [" << subset.lat_start << ":" << lat_end_print << "] = "
    << file_lats[subset.lat_start] << " to " << file_lats[lat_end_print]
    << "  Longitude: [" << subset.lon_start << ":" << lon_end_print << "] = "
    << file_lons[subset.lon_start] << " to " << file_lons[lon_end_print]
    << "  Data reduction: " << subset.lat_count << " x " << subset.lon_count
    << " (was " << nLat << " x " << nLon << ")" << std::endl;
  return subset;
}

// -------------------------------------------------------------------------------------------------
IOStructuredGrid::IOStructuredGrid(const Geometry & geom, const Parameters_ & params)
  : IOBase(geom, params.toConfiguration()), interpolator_(), interpolatorBack_(), geom_(geom),
    gridStr_(""), params_(params), writeFunctionSpace_(), readFunctionSpace_() {
  util::Timer timer(classname(), "IOStructuredGrid");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;

  const std::string mode = params_.mode.value();
  if (!validMode(mode)) {
    ABORT(validModeMessage(mode));
  }

  // Create the Atlas structured grid
  // --------------------------------
  // Create the string to determine the grid name for Atlas
  std::string outputGridType = params.outputGridType.value();

  // Convert the legacy gridtype to what Atlas expects
  if (outputGridType == "latlon") {
    outputGridType = "L" + std::to_string(4*(geom.npx()-1)) + "x" +
                     std::to_string(2*(geom.npy()-1)+1);
  } else if (outputGridType == "gaussian") {
    // Find best matching Gaussian grid
    outputGridType = "F" + std::to_string(geom.npy()-1);
  }

  // Assert that grid begins with either L or F
  if (outputGridType[0] == 'L') {
    gridStr_ = "latlon";
  } else if (outputGridType[0] == 'F') {
    gridStr_ = "gaussian";
  } else {
    // This code is only tested with latlon and regular Gaussian grids. With other grids the code
    // may run but with resulting files containing incorrect or jumbled data.
    ABORT("IOStructuredGrid: outputGridType must begin with L (latlon) or F (regular gaussian). ");
  }

  // Generate the Atlas grid object
  const atlas::Grid grid(outputGridType);

  // Make a custom serial distribution where all points live on rank 0
  // -----------------------------------------------------------------
  std::vector<int> zeros(grid.size(), 0);
  const atlas::grid::Distribution dist(geom.getComm().size(), grid.size(), zeros.data());

  // Create the configuration for the interpolation and populate with the communicator name
  // --------------------------------------------------------------------------------------
  eckit::LocalConfiguration atlas_conf;
  atlas_conf.set("mpi_comm", geom.getComm().name());

  // Structured grid function space
  // ------------------------------
  writeFunctionSpace_.reset(new atlas::functionspace::StructuredColumns(grid, dist, atlas_conf));

  // Create a GeometryData object
  // ----------------------------
  oops::GeometryData geomData(geom.functionSpace(), geom.fields(), geom.levelsAreTopDown(),
                              geom.getComm());

  // -------------------------
  // WRITE-only mode
  // -------------------------
  if (mode == "write") {
    // Create a generic interpolator for converting to the structured grid
    // -------------------------------------------------------------------
    interpolator_.reset(new oops::GlobalInterpolator(params.toConfiguration(),
                                                     geomData,
                                                     *writeFunctionSpace_,
                                                     geom.getComm()));
    oops::Log::trace() << classname() << " constructor done (write)" << std::endl;
    return;
  }

  // -------------------------
  // READ mode: build INPUT geometry from file
  // -------------------------
  if (params_.inputFilename.value() == boost::none) {
    ABORT("IOStructuredGrid: mode is 'read' but no input filename was specified. "
          "The 'input filename' option is required for read mode.");
  }
  const std::string inFile = *params_.inputFilename.value();

  size_t nLat = 0, nLon = 0;
  std::vector<double> file_lats;
  std::vector<double> file_lons;
  if (geom_.getComm().rank() == 0) {
    int ncid;
    nc_rc(nc_open(inFile.c_str(), NC_NOWRITE, &ncid), "nc_open " + inFile);

    const std::string latName = params_.latName.value();
    const std::string lonName = params_.lonName.value();

    int dim_lat, dim_lon;
    nc_rc(nc_inq_dimid(ncid, latName.c_str(), &dim_lat), "nc_inq_dimid " + latName);
    nc_rc(nc_inq_dimlen(ncid, dim_lat, &nLat), "nc_inq_dimlen " + latName);
    nc_rc(nc_inq_dimid(ncid, lonName.c_str(), &dim_lon), "nc_inq_dimid " + lonName);
    nc_rc(nc_inq_dimlen(ncid, dim_lon, &nLon), "nc_inq_dimlen " + lonName);

    file_lats.resize(nLat);
    file_lons.resize(nLon);

    int var_lat, var_lon;
    nc_rc(nc_inq_varid(ncid, latName.c_str(), &var_lat), "nc_inq_varid " + latName);
    nc_rc(nc_get_var_double(ncid, var_lat, file_lats.data()), "nc_get_var_double " + latName);

    nc_rc(nc_inq_varid(ncid, lonName.c_str(), &var_lon), "nc_inq_varid " + lonName);
    nc_rc(nc_get_var_double(ncid, var_lon, file_lons.data()), "nc_get_var_double " + lonName);

    nc_rc(nc_close(ncid), "nc_close");

    oops::Log::info() << "Input file grid: nLon=" << nLon << ", nLat=" << nLat << std::endl;
  }

  // Broadcast dims and coordinates
  geom_.getComm().broadcast(nLat, 0);
  geom_.getComm().broadcast(nLon, 0);
  if (geom_.getComm().rank() != 0) {
    file_lats.resize(nLat);
    file_lons.resize(nLon);
  }
  geom_.getComm().broadcast(file_lats, 0);
  geom_.getComm().broadcast(file_lons, 0);

  // Use an Atlas L-grid only as the structured indexing/distribution container.
  // In HAFS regional applications, the input is treated as a file-defined
  // structured grid using explicit latitude/longitude coordinates from the input file.
  // The exact source coordinates are supplied below through the lonlat field.
  // Native F-grid support can be added later if needed.
  const std::string inputGridType = "L" + std::to_string(nLon) + "x" + std::to_string(nLat);
  const atlas::Grid inputGrid(inputGridType);

  // -----------------------------------------------------------------------------
  // Build distributed structured FunctionSpace for the global file grid
  // -----------------------------------------------------------------------------
  atlas::grid::Partitioner partitioner("equal_regions");

  // --- IMPORTANT: enable halo/ghosts so haloExchange is real and seam stencils exist ---
  atlas_conf.set("halo", 2);

  readFunctionSpace_.reset(
      new atlas::functionspace::StructuredColumns(inputGrid, partitioner, atlas_conf));

  // --- Override lonlat using the configured input longitude/latitude coordinates ---
  atlas::FieldSet structuredAux;

  atlas::Field lonlat = readFunctionSpace_->createField<double>(
      atlas::option::name("lonlat") | atlas::option::variables(2));

  auto lonlatView = atlas::array::make_view<double, 2>(lonlat);
  auto gidxView   = atlas::array::make_view<atlas::gidx_t, 1>(
      readFunctionSpace_->global_index());
  auto ghostView  = atlas::array::make_view<int, 1>(readFunctionSpace_->ghost());

  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (atlas::idx_t p = 0; p < lonlatView.shape(0); ++p) {
    lonlatView(p, 0) = nan;
    lonlatView(p, 1) = nan;
  }

  // The configured latitude coordinate is stored in the same row order as the NetCDF variables.
  // Therefore, the structured row index j is mapped directly to file_lats[j].
  // Do not flip north-to-south files here. Latitude orientation is handled
  // only when computing optional regional subset indices.
  for (atlas::idx_t p = 0; p < lonlatView.shape(0); ++p) {
    if (ghostView(p) != 0) continue;              // owned only
    const auto gidx = gidxView(p);
    if (gidx <= 0) continue;
    const std::size_t g = static_cast<std::size_t>(gidx) - 1;
    const std::size_t i = g % nLon;
    const std::size_t j = g / nLon;
    lonlatView(p, 0) = file_lons[i];
    lonlatView(p, 1) = file_lats[j];
  }

  // populate ghost lonlat consistently
  readFunctionSpace_->haloExchange(lonlat);

  structuredAux.add(lonlat);

  oops::GeometryData structuredGeomData(*readFunctionSpace_,
                                        structuredAux,
                                        geom.levelsAreTopDown(),
                                        geom.getComm());

  // Reverse interpolator (structured -> model)
  eckit::LocalConfiguration interpConfig = params.toConfiguration();
  if (!interpConfig.has("local interpolator type")) {
    interpConfig.set("local interpolator type", "oops unstructured grid interpolator");
  }

  // Only create interpolator once
  if (!interpolatorBack_) {
    interpolatorBack_.reset(new oops::GlobalInterpolator(interpConfig,
                                                         structuredGeomData,
                                                         geom.functionSpace(),
                                                         geom.getComm()));
    if (geom_.getComm().rank() == 0) {
      oops::Log::info() << "Interpolator created (will be reused for all files)" << std::endl;
    }
  }

  oops::Log::trace() << classname() << " constructor done (" << mode << ")" << std::endl;
}
// -------------------------------------------------------------------------------------------------
IOStructuredGrid::~IOStructuredGrid() {
  util::Timer timer(classname(), "~IOStructuredGrid");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}

void IOStructuredGrid::read(State & x,
                            const eckit::LocalConfiguration & fileionames,
                            const eckit::LocalConfiguration & fileioscaling) const {
  util::Timer timer(classname(), "read(State vertical remap)");

  using clock_t = std::chrono::steady_clock;
  auto sec = [](clock_t::time_point a, clock_t::time_point b) {
    return std::chrono::duration<double>(b - a).count();
  };

  const auto t_total0 = clock_t::now();
  const int rank = geom_.getComm().rank();

  auto log0t = [&](const std::string &s, double v) {
    if (rank == 0) oops::Log::info() << s << v << " s" << std::endl;
  };

  const std::string mode = params_.mode.value();
  if (mode == "write") {
    ABORT("IOStructuredGrid::read(State) called with mode='write'. Use mode='read'.");
  }
  if (params_.inputFilename.value() == boost::none) {
    ABORT("IOStructuredGrid::read(State): no input filename was specified. The 'input filename' "
          "option is required for read mode.");
  }

  const std::string inFile = *params_.inputFilename.value();
  const double nan = std::numeric_limits<double>::quiet_NaN();

  // ============================================================
  // 0) Read file dims + coords (rank0) and broadcast
  // ============================================================
  size_t nLat = 0, nLon = 0, nLevFile = 0;
  std::vector<double> file_lats, file_lons;
  int fileNorthToSouth_i = 0;

  auto t_dims0 = clock_t::now();
  if (rank == 0) {
    int ncid;
    nc_rc(nc_open(inFile.c_str(), NC_NOWRITE, &ncid), "nc_open " + inFile);

    const std::string latName = params_.latName.value();
    const std::string lonName = params_.lonName.value();
    const std::string levName = params_.levName.value();

    int dim_lat, dim_lon, dim_lev;
    nc_rc(nc_inq_dimid(ncid, latName.c_str(), &dim_lat), "nc_inq_dimid " + latName);
    nc_rc(nc_inq_dimlen(ncid, dim_lat, &nLat), "nc_inq_dimlen " + latName);
    nc_rc(nc_inq_dimid(ncid, lonName.c_str(), &dim_lon), "nc_inq_dimid " + lonName);
    nc_rc(nc_inq_dimlen(ncid, dim_lon, &nLon), "nc_inq_dimlen " + lonName);

    nc_rc(nc_inq_dimid(ncid, levName.c_str(), &dim_lev), "nc_inq_dimid " + levName);
    nc_rc(nc_inq_dimlen(ncid, dim_lev, &nLevFile), "nc_inq_dimlen " + levName);

    file_lats.resize(nLat);
    file_lons.resize(nLon);

    int var_lat, var_lon;
    nc_rc(nc_inq_varid(ncid, latName.c_str(), &var_lat), "nc_inq_varid " + latName);
    nc_rc(nc_get_var_double(ncid, var_lat, file_lats.data()), "nc_get_var " + latName);
    nc_rc(nc_inq_varid(ncid, lonName.c_str(), &var_lon), "nc_inq_varid " + lonName);
    nc_rc(nc_get_var_double(ncid, var_lon, file_lons.data()), "nc_get_var " + lonName);

    fileNorthToSouth_i = (nLat >= 2 && file_lats[0] > file_lats[nLat - 1]) ? 1 : 0;

    nc_rc(nc_close(ncid), "nc_close " + inFile);
  }

  geom_.getComm().broadcast(nLat, 0);
  geom_.getComm().broadcast(nLon, 0);
  geom_.getComm().broadcast(nLevFile, 0);

  if (rank != 0) {
    file_lats.resize(nLat);
    file_lons.resize(nLon);
  }
  geom_.getComm().broadcast(file_lats.begin(), file_lats.end(), 0);
  geom_.getComm().broadcast(file_lons.begin(), file_lons.end(), 0);
  geom_.getComm().broadcast(fileNorthToSouth_i, 0);

  log0t("[TIMER] read dims+coords+bcast: ", sec(t_dims0, clock_t::now()));

  // ============================================================
  // 1) Optional regional subset (lon/lat bounds) to reduce I/O
  // ============================================================
  size_t lat_start = 0, lat_count = nLat, lon_start = 0, lon_count = nLon;
  const size_t halo = 10;

  int use_regional_subset_i = 0;
  double lon_min = -180.0, lon_max = 180.0, lat_min = -90.0, lat_max = 90.0;

  if (params_.lon_min.value() != boost::none) {
    lon_min = *params_.lon_min.value();
    use_regional_subset_i = 1;
  }
  if (params_.lon_max.value() != boost::none) {
    lon_max = *params_.lon_max.value();
    use_regional_subset_i = 1;
  }
  if (params_.lat_min.value() != boost::none) {
    lat_min = *params_.lat_min.value();
    use_regional_subset_i = 1;
  }
  if (params_.lat_max.value() != boost::none) {
    lat_max = *params_.lat_max.value();
    use_regional_subset_i = 1;
  }

  const bool fileNorthToSouth = (fileNorthToSouth_i != 0);
  const double lat_lo_user = std::min(lat_min, lat_max);
  const double lat_hi_user = std::max(lat_min, lat_max);
  const double lon_lo_user = std::min(lon_min, lon_max);
  const double lon_hi_user = std::max(lon_min, lon_max);

  auto t_subset0 = clock_t::now();
  if (rank == 0 && use_regional_subset_i == 1) {
    const RegionalSubset subset = computeRegionalSubset(
        file_lats, file_lons, nLat, nLon, fileNorthToSouth,
        lat_lo_user, lat_hi_user, lon_lo_user, lon_hi_user, halo);
    lat_start = subset.lat_start;
    lat_count = subset.lat_count;
    lon_start = subset.lon_start;
    lon_count = subset.lon_count;
  }
  log0t("[TIMER] subset compute: ", sec(t_subset0, clock_t::now()));

  geom_.getComm().broadcast(use_regional_subset_i, 0);
  const bool use_regional_subset = (use_regional_subset_i != 0);
  geom_.getComm().broadcast(lat_start, 0);
  geom_.getComm().broadcast(lat_count, 0);
  geom_.getComm().broadcast(lon_start, 0);
  geom_.getComm().broadcast(lon_count, 0);

  const size_t lat_end = lat_start + lat_count - 1;
  const size_t lon_end = lon_start + lon_count - 1;

  // ============================================================
  // 2) Pull State FieldSet and get model levels
  // ============================================================
  atlas::FieldSet fieldsModelAll;
  x.toFieldSet(fieldsModelAll);

  if (!fieldsModelAll.has("air_temperature")) {
    ABORT("State missing air_temperature (needed to infer nLevModel)");
  }
  atlas::Field & tField = fieldsModelAll.field("air_temperature");
  const int nLevModel = (tField.rank() >= 2) ? static_cast<int>(tField.shape(1)) : 1;

  if (rank == 0) {
    oops::Log::info() << "[CHECK-NLEV] file " << params_.levName.value() << "=" << nLevFile
                      << " model levels=" << nLevModel << std::endl;
  }
  if (nLevModel <= 1) ABORT("Model State appears not to have 3D levels (nLevModel<=1)");

  // ============================================================
  // 3) Get target ak/bk from the FV3 Geometry
  // ============================================================
  const std::vector<double> & ak = geom_.ak();
  const std::vector<double> & bk = geom_.bk();
  const size_t nLevModelP1 = static_cast<size_t>(nLevModel + 1);
  if (ak.size() != nLevModelP1 || bk.size() != nLevModelP1) {
    if (rank == 0) {
      oops::Log::error() << "Geometry ak/bk size mismatch: expected "
                         << nLevModelP1 << " got ak=" << ak.size()
                         << " bk=" << bk.size() << std::endl;
    }
    ABORT("Geometry ak/bk size mismatch");
  }

  // ============================================================
  // 4) Variable mapping
  // ============================================================
  struct Var3D { std::string stateName; std::string fileName; };
  std::vector<Var3D> vars3d;

  const auto & stateVars = x.variables();

  std::map<std::string, std::string> ioNameMap;
  if (fileionames.has("field io names")) {
    eckit::LocalConfiguration ionames = fileionames.getSubConfiguration("field io names");
    std::vector<std::string> keys = ionames.keys();
    for (const auto & key : keys) {
      ioNameMap[key] = ionames.getString(key);
    }
  } else {
    std::vector<std::string> keys = fileionames.keys();
    for (const auto & key : keys) {
      try {
        ioNameMap[key] = fileionames.getString(key);
      } catch (...) {
      }
    }
  }

  const std::vector<std::string> skip2D = {"air_pressure_at_surface", "surface_pressure",
                                      "geopotential_height_at_surface"};
  for (size_t i = 0; i < stateVars.size(); ++i) {
    std::string varName = stateVars[i].name();
    if (std::find(skip2D.begin(), skip2D.end(), varName) != skip2D.end()) continue;
    std::string fileName = varName;
    if (ioNameMap.count(varName) > 0) fileName = ioNameMap[varName];
    vars3d.push_back({varName, fileName});
  }

  if (rank == 0) {
    oops::Log::info() << "[I/O] Reading " << vars3d.size() << " 3D variables:" << std::endl;
    for (const auto & v : vars3d) {
      oops::Log::info() << "  " << v.stateName << " <- " << v.fileName << std::endl;
    }
  }

  for (const auto & v : vars3d) {
    if (!fieldsModelAll.has(v.stateName)) {
      if (rank == 0) {
        oops::Log::error() << "State missing required field: " << v.stateName << std::endl;
      }
      ABORT("State missing required 3D field");
    }
  }
  if (!fieldsModelAll.has("air_pressure_at_surface")) {
    ABORT("State missing air_pressure_at_surface");
  }

  const std::string dpresName = "__pressure_thickness__";
  const std::string psName    = "__surface_pressure__";

  auto requiredFileName = [&](const std::string & fieldName,
                              const std::string & purpose) -> std::string {
    const auto it = ioNameMap.find(fieldName);
    if (it == ioNameMap.end()) {
      ABORT("IOStructuredGrid::read(State): field io names must define '" + fieldName
            + "' for " + purpose + ".");
    }
    return it->second;
  };

  const std::string dpresFileName =
      requiredFileName("air_pressure_thickness", "source pressure thickness");
  const std::string psFileName =
      requiredFileName("air_pressure_at_surface", "source surface pressure");

  const std::string orogName = "geopotential_height_at_surface";
  const bool doVertRemap = params_.doVerticalRemapping.value();
  std::string sourceOrogFileName;
  atlas::FieldSet fieldsOrog;

  if (doVertRemap) {
    sourceOrogFileName = requiredFileName(orogName, "source orography for VertRemap");

    int ncid_check;
    nc_rc(nc_open(inFile.c_str(), NC_NOWRITE, &ncid_check),
          "nc_open for source orography check " + inFile);
    int varid_check;
    const int orog_rc = nc_inq_varid(ncid_check, sourceOrogFileName.c_str(), &varid_check);
    nc_rc(nc_close(ncid_check), "nc_close source orography check");
    if (orog_rc != NC_NOERR) {
      ABORT("IOStructuredGrid::read(State): source orography variable '"
            + sourceOrogFileName + "' was requested for VertRemap but is not present in "
            + inFile + ".");
    }

    if (geom_.fields().has(orogName)) {
      atlas::Field zsOrogTarget = geom_.fields().field(orogName).clone();
      zsOrogTarget.rename(orogName);
      fieldsOrog.add(zsOrogTarget);
    } else if (fieldsModelAll.has(orogName)) {
      atlas::Field zsOrogTarget = fieldsModelAll.field(orogName).clone();
      zsOrogTarget.rename(orogName);
      fieldsOrog.add(zsOrogTarget);
    } else {
      ABORT("IOStructuredGrid::read(State): do vertical remapping is true, but target field '"
            + orogName + "' is missing from Geometry fields and State FieldSet.");
    }

    if (rank == 0) {
      oops::Log::info() << "[VERT-REMAP] VertRemap enabled in read mode using source "
                        << sourceOrogFileName << " and target " << orogName << std::endl;
    }
  } else if (rank == 0) {
    oops::Log::info() << "[VERT-REMAP] VertRemap disabled in read mode; "
                      << "using pressure-space log-p remap only." << std::endl;
  }

  // ============================================================
  // 5) Create source fields on structured FS and temporary model-grid
  //    fields still at FILE vertical levels. Vertical remap will be
  //    done AFTER horizontal interpolation.
  // ============================================================
  atlas::FieldSet srcFS;
  atlas::FieldSet tgtModelFileLevels;

  auto make_model_tmp = [&](const std::string & name, int levels) -> atlas::Field {
    return geom_.functionSpace().createField<double>(
      atlas::option::name(name) | atlas::option::levels(levels));
  };

  for (const auto & v : vars3d) {
    srcFS.add(readFunctionSpace_->createField<double>(
      atlas::option::name(v.stateName) | atlas::option::levels(static_cast<int>(nLevFile))));
    tgtModelFileLevels.add(make_model_tmp(v.stateName, static_cast<int>(nLevFile)));
  }

  srcFS.add(readFunctionSpace_->createField<double>(
    atlas::option::name(dpresName) | atlas::option::levels(static_cast<int>(nLevFile))));
  tgtModelFileLevels.add(make_model_tmp(dpresName, static_cast<int>(nLevFile)));

  srcFS.add(readFunctionSpace_->createField<double>(
    atlas::option::name(psName) | atlas::option::levels(1)));
  tgtModelFileLevels.add(make_model_tmp(psName, 1));

  if (doVertRemap) {
    srcFS.add(readFunctionSpace_->createField<double>(
      atlas::option::name(orogName) | atlas::option::levels(1)));
    tgtModelFileLevels.add(make_model_tmp(orogName, 1));
  }

  auto set_interp_type = [&](atlas::FieldSet & fs, const std::string & val) {
    for (auto & f : fs) {
      if (!f.metadata().has("interp_type")) f.metadata().set("interp_type", val);
    }
  };
  set_interp_type(srcFS, "default");
  set_interp_type(tgtModelFileLevels, "default");

  // ============================================================
  // 6) Read slabs (level-parallel I/O) into structured srcFS
  // ============================================================
  const auto gidxView = atlas::array::make_view<atlas::gidx_t, 1>(
      readFunctionSpace_->global_index());
  const auto ghostView = atlas::array::make_view<int, 1>(readFunctionSpace_->ghost());

  const int n_ranks = geom_.getComm().size();

  struct VarInfo {
    std::string fileName;
    std::string stateName;
    int nLevels;
    bool is3D;
  };

  std::vector<VarInfo> all_vars;
  for (const auto & v : vars3d) {
    all_vars.push_back({v.fileName, v.stateName, static_cast<int>(nLevFile), true});
  }
  all_vars.push_back({dpresFileName, dpresName, static_cast<int>(nLevFile), true});
  all_vars.push_back({psFileName, psName, 1, false});
  if (doVertRemap) {
    all_vars.push_back({sourceOrogFileName, orogName, 1, false});
  }

  struct SliceWork {
    int var_idx;
    int level_idx;
  };

  std::vector<SliceWork> all_slices;
  for (int v = 0; v < static_cast<int>(all_vars.size()); ++v) {
    for (int lev = 0; lev < all_vars[v].nLevels; ++lev) {
      all_slices.push_back({v, lev});
    }
  }

  const int total_slices = static_cast<int>(all_slices.size());
  const int slices_per_rank = (total_slices + n_ranks - 1) / n_ranks;

  if (rank == 0) {
    oops::Log::info() << "[I/O LEVEL-PARALLEL] " << total_slices
                      << " 2D slices distributed across " << n_ranks << " ranks ("
                      << slices_per_rank << " slices/rank avg)" << std::endl;
  }

  int ncid;
  nc_rc(nc_open(inFile.c_str(), NC_NOWRITE, &ncid), "nc_open " + inFile);

  auto t_read_start = clock_t::now();
  std::map<int, std::map<int, std::vector<float>>> my_slices;
  const size_t plane_size = lat_count * lon_count;

  int my_slice_count = 0;
  for (int s = rank; s < total_slices; s += n_ranks) {
    const auto & work = all_slices[s];
    const int v_idx = work.var_idx;
    const int lev_idx = work.level_idx;
    const auto & var = all_vars[v_idx];

    std::vector<float> slice_data(plane_size);
    int varid;
    nc_rc(nc_inq_varid(ncid, var.fileName.c_str(), &varid), "nc_inq_varid " + var.fileName);

    if (var.is3D) {
      size_t start[4] = {0, static_cast<size_t>(lev_idx), lat_start, lon_start};
      size_t count[4] = {1, 1, lat_count, lon_count};
      nc_rc(nc_get_vara_float(ncid, varid, start, count, slice_data.data()),
            "nc_get_vara_float " + var.fileName);
    } else {
      size_t start[3] = {0, lat_start, lon_start};
      size_t count[3] = {1, lat_count, lon_count};
      nc_rc(nc_get_vara_float(ncid, varid, start, count, slice_data.data()),
            "nc_get_vara_float " + var.fileName);
    }

    my_slices[v_idx][lev_idx] = std::move(slice_data);
    my_slice_count++;
  }

  nc_rc(nc_close(ncid), "nc_close");
  double t_read_local = sec(t_read_start, clock_t::now());

  if (rank == 0 || rank == n_ranks - 1) {
    oops::Log::info() << "[I/O rank " << rank << "] read " << my_slice_count
                      << " slices in " << t_read_local << "s" << std::endl;
  }

  auto t_gather_start = clock_t::now();
  for (int v = 0; v < static_cast<int>(all_vars.size()); ++v) {
    const auto & var = all_vars[v];
    std::vector<float> full_var(plane_size * var.nLevels);

    int var_slice_offset = 0;
    for (int vv = 0; vv < v; ++vv) var_slice_offset += all_vars[vv].nLevels;

    for (int lev = 0; lev < var.nLevels; ++lev) {
      const int slice_global_idx = var_slice_offset + lev;
      const int owner_rank = slice_global_idx % n_ranks;
      std::vector<float> level_data(plane_size);
      if (rank == owner_rank) level_data = my_slices[v][lev];
      geom_.getComm().broadcast(level_data.begin(), level_data.end(), owner_rank);
      std::copy(level_data.begin(), level_data.end(), full_var.begin() + lev * plane_size);
    }

    atlas::Field & field = srcFS.field(var.stateName);
    auto v_view = atlas::array::make_view<double, 2>(field);
    const int nLev = var.nLevels;

    for (atlas::idx_t p = 0; p < v_view.shape(0); ++p) {
      for (int k = 0; k < nLev; ++k) v_view(p, k) = nan;
    }

    for (atlas::idx_t p = 0; p < v_view.shape(0); ++p) {
      if (ghostView(p) != 0) continue;
      const auto gidx = gidxView(p);
      if (gidx <= 0) continue;

      const std::size_t g = static_cast<std::size_t>(gidx) - 1;
      const std::size_t i = g % nLon;
      const std::size_t j = g / nLon;
      if (i >= nLon || j >= nLat) continue;

      // The global index row j is used as the raw file row. This must remain
      // consistent with the lonlat assignment above: row j uses file_lats[j]
      // and data row j.
      if (use_regional_subset) {
        if (i < lon_start || i > lon_end || j < lat_start || j > lat_end) continue;
      }

      const std::size_t ii = i - lon_start;
      const std::size_t jj = j - lat_start;
      for (int k = 0; k < nLev; ++k) {
        const size_t idx = (static_cast<size_t>(k) * lat_count + jj) * lon_count + ii;
        v_view(p, k) = static_cast<double>(full_var[idx]);
      }
    }

    readFunctionSpace_->haloExchange(field);
  }

  double t_gather = sec(t_gather_start, clock_t::now());
  double t_read_max = t_read_local;
  geom_.getComm().allReduceInPlace(t_read_max, eckit::mpi::Operation::MAX);
  log0t("[TIMER LEVEL-PARALLEL] NetCDF read (parallel across levels, max): ", t_read_max);
  log0t("[TIMER LEVEL-PARALLEL] Gather and assemble (broadcasts): ", t_gather);

  // ============================================================
  // 7) Horizontal interpolation first: file grid -> model grid,
  //    still preserving FILE vertical levels.
  // ============================================================
  auto t_h0 = clock_t::now();

  for (auto & f : tgtModelFileLevels) {
    auto v = atlas::array::make_view<double, 2>(f);
    for (atlas::idx_t p = 0; p < v.shape(0); ++p) {
      for (int k = 0; k < v.shape(1); ++k) v(p, k) = nan;
    }
  }

  interpolatorBack_->apply(srcFS, tgtModelFileLevels);
  log0t("[TIMER] horizontal interp (file->model @file-levels): ", sec(t_h0, clock_t::now()));

  // ============================================================
  // 8) Vertical remap on MODEL grid: file levels -> model levels
  // ============================================================
  auto t_v0 = clock_t::now();

  atlas::FieldSet tgtModelRemapped;
  for (const auto & v : vars3d) {
    tgtModelRemapped.add(make_model_tmp(v.stateName, nLevModel));
  }
  set_interp_type(tgtModelRemapped, "default");

  auto dpresV_model_fileLev =
      atlas::array::make_view<double, 2>(tgtModelFileLevels.field(dpresName));
  auto psV_model = atlas::array::make_view<double, 2>(tgtModelFileLevels.field(psName));
  auto psState =
      atlas::array::make_view<double, 2>(fieldsModelAll.field("air_pressure_at_surface"));

  auto interp_logp = [&](const std::vector<double> & p_src,
                         const std::vector<double> & x_src,
                         double p_tgt) -> double {
    const size_t n = p_src.size();
    if (n < 2) return x_src.empty() ? nan : x_src.front();
    if (p_tgt <= p_src.front()) return x_src.front();
    if (p_tgt >= p_src.back())  return x_src.back();

    auto it = std::upper_bound(p_src.begin(), p_src.end(), p_tgt);
    size_t k1 = std::max<size_t>(1, static_cast<size_t>(it - p_src.begin())) - 1;
    size_t k2 = k1 + 1;

    const double p1 = std::max(1.0, p_src[k1]);
    const double p2 = std::max(1.0, p_src[k2]);
    const double x1 = x_src[k1];
    const double x2 = x_src[k2];
    const double denom = std::log(p2) - std::log(p1);
    if (std::abs(denom) < 1.0e-12) return x1;

    const double w = (std::log(std::max(1.0, p_tgt)) - std::log(p1)) / denom;
    return x1 + w * (x2 - x1);
  };

  std::vector<double> p_int_src(nLevFile + 1);
  std::vector<double> p_mid_src(nLevFile);
  std::vector<double> p_int_tgt(nLevModel + 1);
  std::vector<double> p_mid_tgt(nLevModel);
  std::vector<double> x_src(nLevFile);

  const atlas::idx_t npts_model = psV_model.shape(0);
  for (const auto & vinfo : vars3d) {
    auto srcVar_fileLev =
        atlas::array::make_view<double, 2>(tgtModelFileLevels.field(vinfo.stateName));
    auto dstVar_modLev =
        atlas::array::make_view<double, 2>(tgtModelRemapped.field(vinfo.stateName));

    for (atlas::idx_t p = 0; p < npts_model; ++p) {
      p_int_src[0] = 0.0;
      for (size_t k = 0; k < nLevFile; ++k) {
        double dp = dpresV_model_fileLev(p, static_cast<int>(k));
        if (!std::isfinite(dp) || dp < 0.0) dp = 0.0;
        p_int_src[k + 1] = p_int_src[k] + dp;
      }

      double ps_col = psV_model(p, 0);
      if ((!std::isfinite(ps_col) || ps_col <= 0.0) && p_int_src[nLevFile] > 0.0) {
        ps_col = p_int_src[nLevFile];
      }
      if ((!std::isfinite(ps_col) || ps_col <= 0.0) &&
          std::isfinite(psState(p, 0)) && psState(p, 0) > 0.0) {
        ps_col = psState(p, 0);
      }

      const double sumdp = p_int_src[nLevFile];
      const double scale = (sumdp > 0.0 && ps_col > 0.0) ? (ps_col / sumdp) : 1.0;
      for (size_t k = 0; k <= nLevFile; ++k) p_int_src[k] *= scale;

      for (size_t k = 0; k < nLevFile; ++k) {
        p_mid_src[k] = 0.5 * (p_int_src[k] + p_int_src[k + 1]);
        if (p_mid_src[k] < 1.0) p_mid_src[k] = 1.0;
        x_src[k] = srcVar_fileLev(p, static_cast<int>(k));
      }

      for (int k = 0; k <= nLevModel; ++k) {
        p_int_tgt[k] = ak[k] + bk[k] * ps_col;
        if (p_int_tgt[k] < 1.0) p_int_tgt[k] = 1.0;
      }
      for (int k = 0; k < nLevModel; ++k) {
        p_mid_tgt[k] = 0.5 * (p_int_tgt[k] + p_int_tgt[k + 1]);
        if (p_mid_tgt[k] < 1.0) p_mid_tgt[k] = 1.0;
      }

      for (int k = 0; k < nLevModel; ++k) {
        dstVar_modLev(p, k) = interp_logp(p_mid_src, x_src, p_mid_tgt[k]);
      }
    }
  }

  if (doVertRemap) {
    // Match the write-side usage pattern:
    //   - tgtModelRemapped carries the horizontally interpolated source surface height
    //     as geopotential_height_at_surface.
    //   - fieldsOrog carries the real target-grid/model surface height.
    //   - VertRemap performs its own terrain-aware surface-pressure adjustment and
    //     vertical remapping. No local simplified hydrostatic formula is used here.
    if (!tgtModelRemapped.has(orogName)) {
      atlas::Field zsSourceInterp = tgtModelFileLevels.field(orogName).clone();
      zsSourceInterp.rename(orogName);
      tgtModelRemapped.add(zsSourceInterp);
    }
    if (!tgtModelRemapped.has("air_pressure_at_surface")) {
      // This is the horizontally interpolated source surface pressure. It is
      // intentionally added to the FieldSet passed to VertRemap. VertRemap may
      // replace it with terrain-adjusted surface pressure, and that adjusted
      // value is propagated back to the State below.
      atlas::Field psSourceInterp = make_model_tmp("air_pressure_at_surface", 1);
      auto psSrc = atlas::array::make_view<double, 2>(tgtModelFileLevels.field(psName));
      auto psDst = atlas::array::make_view<double, 2>(psSourceInterp);
      for (atlas::idx_t p = 0; p < psDst.shape(0); ++p) {
        psDst(p, 0) = psSrc(p, 0);
      }
      tgtModelRemapped.add(psSourceInterp);
    }

    fv3jedi::VertRemap vert_remap(geom_, fieldsOrog);
    tgtModelRemapped = vert_remap.remap(tgtModelRemapped);
  }

  log0t("[TIMER] vertical remap (model-grid, file-levels->model-levels): ",
        sec(t_v0, clock_t::now()));

  // ============================================================
  // 9) Copy remapped fields into State and update surface pressure
  // ============================================================
  for (const auto & vinfo : vars3d) {
    auto src = atlas::array::make_view<double, 2>(tgtModelRemapped.field(vinfo.stateName));
    auto dst = atlas::array::make_view<double, 2>(fieldsModelAll.field(vinfo.stateName));
    for (atlas::idx_t p = 0; p < dst.shape(0); ++p) {
      for (int k = 0; k < nLevModel; ++k) dst(p, k) = src(p, k);
    }
  }

  if (tgtModelRemapped.has("air_pressure_at_surface")) {
    // In the VertRemap path this is the terrain-adjusted surface pressure
    // returned by VertRemap. If VertRemap is disabled, this field is absent, so
    // the State receives the horizontally interpolated source surface pressure below.
    auto psRemapped = atlas::array::make_view<double, 2>(
      tgtModelRemapped.field("air_pressure_at_surface"));
    for (atlas::idx_t p = 0; p < psState.shape(0); ++p) {
      const double v = psRemapped(p, 0);
      if (std::isfinite(v) && v > 0.0) psState(p, 0) = v;
    }
  } else {
    for (atlas::idx_t p = 0; p < psState.shape(0); ++p) {
      const double v = psV_model(p, 0);
      if (std::isfinite(v) && v > 0.0) psState(p, 0) = v;
    }
  }

  // ============================================================
  // 10) Copy back to State
  // ============================================================
  auto t_from0 = clock_t::now();
  x.fromFieldSet(fieldsModelAll);
  log0t("[TIMER] fromFieldSet: ", sec(t_from0, clock_t::now()));

  log0t("[TIMER] TOTAL read(): ", sec(t_total0, clock_t::now()));
}
// -------------------------------------------------------------------------------------------------
void IOStructuredGrid::read(Increment & dx, const eckit::LocalConfiguration & fileionames,
                                const eckit::LocalConfiguration & fileioscaling) const {
  ABORT("IOStructuredGrid::read(Increment) not implemented");
}

// -------------------------------------------------------------------------------------------------

template <typename T>
void IOStructuredGrid::interpAndWrite(const T & obj, const std::string & label,
                                      const eckit::LocalConfiguration & fileionames,
                                      const eckit::LocalConfiguration & fileioscaling) const {
  util::Timer timer(classname(), "write " + label);
  oops::Log::trace() << classname() << " write " << label << " starting" << std::endl;

  // Create field sets
  atlas::FieldSet fieldsCubeSphere;
  atlas::FieldSet fieldsGeographic;
  obj.toFieldSet(fieldsCubeSphere);

  // Apply interpolation
  interpolator_->apply(fieldsCubeSphere, fieldsGeographic);

  // If orography filename is provided, remap vertical coordinates
  if ( !params_.doVerticalRemapping.value() ) {
    // Write to disk if rank 0
    if (geom_.getComm().rank() == 0) {
      this->writeStructuredFields(fieldsGeographic, obj.validTime(), fileionames, fileioscaling);
    }
  } else {
    ASSERT(params_.orographyFilename.value() != boost::none);

    // Define orography variables
    atlas::FieldSet fieldsOrog;
    atlas::Field zsOrogNew = fieldsGeographic["geopotential_height_at_surface"].clone();
    fieldsOrog.add(zsOrogNew);

    // Write to disk if rank 0
    if (geom_.getComm().rank() == 0) {
      // Read structured-grid orography from file
      const std::string orogFilename = params_.orographyFilename.value().value();
      this -> readStructuredFields(orogFilename, fieldsOrog, obj.validTime(),
                                   fileionames, fileioscaling);

      // Remap the vertical coordinates to account for orography changes
      fv3jedi::VertRemap vert_remap(geom_, fieldsOrog);
      atlas::FieldSet fieldsGeographicRemap = vert_remap.remap(fieldsGeographic);

      // Write to disk
      this->writeStructuredFields(fieldsGeographicRemap, obj.validTime(),
                                  fileionames, fileioscaling);
    }
  }

  oops::Log::trace() << classname() << " write " << label << " done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void IOStructuredGrid::write(const State & x, const eckit::LocalConfiguration & fileionames,
                                const eckit::LocalConfiguration & fileioscaling) const {
  this->interpAndWrite(x, "state", fileionames, fileioscaling);
}

// -------------------------------------------------------------------------------------------------

void IOStructuredGrid::write(const Increment & dx, const eckit::LocalConfiguration & fileionames,
                             const eckit::LocalConfiguration & fileioscaling) const {
  this->interpAndWrite(dx, "increment", fileionames, fileioscaling);
}

// -------------------------------------------------------------------------------------------------

void IOStructuredGrid::writeStructuredFields(const atlas::FieldSet & fields,
                                             const util::DateTime & time,
                                             const eckit::LocalConfiguration & ioNames,
                                             const eckit::LocalConfiguration & ioScaling) const {
  // NetCDF IDs
  // ----------
  int fileId;

  // Dimension indices
  int latId;
  int lonId;
  int levId;
  int edgId;
  int forId;
  int timId;

  // Variable indices
  int fIv;
  std::map<std::string, int> fieldIvs;

  // Get ak/bk for writing
  // ---------------------
  std::vector<double> ak = geom_.ak();
  std::vector<double> bk = geom_.bk();

  // Get the name of the file and adjust with datetime
  // -------------------------------------------------
  std::string pathFile = params_.filename.value();

  // For backward compatibility add some things to the filename if not already present
  if (pathFile.find("%Y") == std::string::npos) {
    pathFile += "%Y%m%d_%H%M%Sz";
  }
  if (pathFile.find(".nc") == std::string::npos) {
    pathFile += ".nc4";
  }

  // Format the datetime string
  pathFile = time.formatString(pathFile);

  // Replace member number (ensemble applciaitons)
  util::stringfunctions::swapNameMember(params_.toConfiguration(), pathFile);

  // Create a file to write fields into
  // ----------------------------------
  nc_rc(nc_create(pathFile.c_str(), NC_CLOBBER | NC_NETCDF4, &fileId), "nc_create" + pathFile);

  // Create regular grid for determining lat/lon values
  // --------------------------------------------------
  const atlas::RegularGrid regGrid(writeFunctionSpace_->grid());

  // Define the dimensions in the file
  // ---------------------------------
  const int nLat = regGrid.ny();
  const int nLon = regGrid.nx();
  const int nLev = geom_.npz();
  const int nEdg = geom_.npz() + 1;
  const int nFor = 4;
  const int nTim = 1;

  nc_rc(nc_def_dim(fileId, params_.latName.value().c_str(), nLat, &latId), "nc_def_dim (lat)");
  nc_rc(nc_def_dim(fileId, params_.lonName.value().c_str(), nLon, &lonId), "nc_def_dim (lon)");
  nc_rc(nc_def_dim(fileId, params_.levName.value().c_str(), nLev, &levId), "nc_def_dim (lev)");
  nc_rc(nc_def_dim(fileId, params_.edgName.value().c_str(), nEdg, &edgId), "nc_def_dim (edg)");
  nc_rc(nc_def_dim(fileId, params_.forName.value().c_str(), nFor, &forId), "nc_def_dim (for)");
  nc_rc(nc_def_dim(fileId, params_.timName.value().c_str(), nTim, &timId), "nc_def_dim (tim)");

  // Define the dimensions variables in the file
  // -------------------------------------------
  std::vector<double> latArr(nLat);
  std::vector<double> lonArr(nLon);
  std::vector<int> levArr(nLev);
  std::vector<int> edgArr(nEdg);
  std::vector<int> forArr(nFor);
  std::vector<int> timArr(nTim);

  for (int i = 0; i < nLat; ++i) {
    latArr[i] = regGrid.y(nLat - 1 - i);
  }
  for (int i = 0; i < nLon; ++i) {
    lonArr[i] = regGrid.x(i);
  }
  for (int i = 0; i < nLev; ++i) {
    levArr[i] = i + 1;
  }
  for (int i = 0; i < nEdg; ++i) {
    edgArr[i] = i + 1;
  }
  for (int i = 0; i < nFor; ++i) {
    forArr[i] = i + 1;
  }
  for (int i = 0; i < nTim; ++i) {
    timArr[i] = i + 1;
  }

  // Write the dimension variables (and attributes) to the file
  // ----------------------------------------------------------
  nc_rc(nc_def_var(fileId, params_.latName.value().c_str(), NC_DOUBLE, 1, &latId, &fIv),
        "nc_def_var (lat)");
  nc_rc(nc_put_att_text(fileId, fIv, "units", strlen("degrees_north"), "degrees_north"),
        "nc_put_att_text (lat)");
  fieldIvs[params_.latName.value()] = fIv;

  nc_rc(nc_def_var(fileId, params_.lonName.value().c_str(), NC_DOUBLE, 1, &lonId, &fIv),
        "nc_def_var (lon)");
  nc_rc(nc_put_att_text(fileId, fIv, "units", strlen("degrees_east"), "degrees_east"),
        "nc_put_att_text (lon)");
  fieldIvs[params_.lonName.value()] = fIv;

  nc_rc(nc_def_var(fileId, params_.levName.value().c_str(), NC_INT, 1, &levId, &fIv),
        "nc_def_var (lev)");
  nc_rc(nc_put_att_text(fileId, fIv, "units", strlen("1"), "1"),
        "nc_put_att_text (lev)");
  fieldIvs[params_.levName.value()] = fIv;

  nc_rc(nc_def_var(fileId, params_.edgName.value().c_str(), NC_INT, 1, &edgId, &fIv),
        "nc_def_var (edg)");
  nc_rc(nc_put_att_text(fileId, fIv, "units", strlen("1"), "1"),
        "nc_put_att_text (edg)");
  fieldIvs[params_.edgName.value()] = fIv;

  nc_rc(nc_def_var(fileId, params_.forName.value().c_str(), NC_INT, 1, &forId, &fIv),
        "nc_def_var (for)");
  nc_rc(nc_put_att_text(fileId, fIv, "units", strlen("1"), "1"),
        "nc_put_att_text (for)");
  fieldIvs[params_.forName.value()] = fIv;

  nc_rc(nc_def_var(fileId, params_.timName.value().c_str(), NC_INT, 1, &timId, &fIv),
        "nc_def_var (tim)");
  nc_rc(nc_put_att_text(fileId, fIv, "units", strlen("1"), "1"),
        "nc_put_att_text (tim)");
  fieldIvs[params_.timName.value()] = fIv;

  // Define some categories of dimension IDs for fields
  // --------------------------------------------------
  std::map<int, std::vector<int>> fieldDims;
  fieldDims[nLev] = {timId, levId, latId, lonId};  // Fields at levels
  fieldDims[nEdg] = {timId, edgId, latId, lonId};  // Fields at edges
  fieldDims[4] = {timId, forId, latId, lonId};     // Fields at four levels
  fieldDims[1] = {timId, latId, lonId};            // Fields at surface
  fieldDims[0] = {timId, latId, lonId};            // Fields at surface

  // Set float precision for fields
  // ------------------------------
  const int floatPrecision = params_.floatPrecision.value();
  const int ncPrec = (floatPrecision == 4) ? NC_FLOAT : NC_DOUBLE;

  // Define all the fields that will be written
  // ------------------------------------------
  for (auto& field : fields) {
    // Get number of levels for this field
    const int nLevField = field.shape(1);

    // Get dimensions for this field from map
    auto it = fieldDims.find(nLevField);
    if (it == fieldDims.end()) {
      std::ostringstream oss;
      oss << "IOStructuredGrid::writeStructuredFields: "
          << "No entry in fieldDims for field '" << field.name()
          << "' with " << nLevField << " levels.";
      ABORT(oss.str());
    }
    const auto &dims = it->second;

    // Look for fieldname in the iofile configuration and use the value if key found
    const std::string fieldLong = field.name();
    const char * fieldLongC = fieldLong.c_str();

    // Get the fieldmetadata for this field
    const FieldMetadata & fieldMetadata = geom_.fieldsMetaData().getFieldMetadata(fieldLong);
    std::string unitsStr = fieldMetadata.getVarUnits();
    const char * units = unitsStr.c_str();

    std::string fieldName = fieldLong;
    if (ioNames.has(fieldName)) {
      fieldName = ioNames.getString(fieldLong);
    }

    // Define the field in the file
    nc_rc(nc_def_var(fileId, fieldName.c_str(), ncPrec, dims.size(), dims.data(), &fIv),
          "nc_def_var " + fieldName);
    nc_rc(nc_put_att_text(fileId, fIv, "units", strlen(units), units),
          "nc_put_att_text " + fieldName + " units");
    nc_rc(nc_put_att_text(fileId, fIv, "long_name", strlen(fieldLongC), fieldLongC),
          "nc_put_att_text " + fieldName + " long_name");

    // Insert field into the fieldIvs map
    fieldIvs[field.name()] = fIv;
  }

  // Write ak/bk to the file as global attributes
  // --------------------------------------------
  nc_rc(nc_put_att_double(fileId, NC_GLOBAL, "ak", NC_DOUBLE, ak.size(), ak.data()),
          "nc_put_att_double (ak)");
  nc_rc(nc_put_att_double(fileId, NC_GLOBAL, "bk", NC_DOUBLE, bk.size(), bk.data()),
          "nc_put_att_double (bk)");
  nc_rc(nc_put_att_text(fileId, NC_GLOBAL, "grid", strlen(gridStr_.c_str()), gridStr_.c_str()),
          "nc_put_att_text (grid)");
  nc_rc(nc_put_att_int(fileId, NC_GLOBAL, "im", NC_INT, 1, &nLon),
          "nc_put_att_int (im)");
  nc_rc(nc_put_att_int(fileId, NC_GLOBAL, "jm", NC_INT, 1, &nLat),
          "nc_put_att_int (im)");

  // End definition mode
  // -------------------
  nc_rc(nc_enddef(fileId), "nc_enddef");

  // Write coordinate data into the file
  // -----------------------------------
  nc_rc(nc_put_var_double(fileId, fieldIvs[params_.latName.value()], latArr.data()),
        "nc_put_var_double (lat)");
  nc_rc(nc_put_var_double(fileId, fieldIvs[params_.lonName.value()], lonArr.data()),
        "nc_put_var_double (lon)");
  nc_rc(nc_put_var_int(fileId, fieldIvs[params_.levName.value()], levArr.data()),
        "nc_put_var_int (lev)");
  nc_rc(nc_put_var_int(fileId, fieldIvs[params_.edgName.value()], edgArr.data()),
        "nc_put_var_int (edg)");
  nc_rc(nc_put_var_int(fileId, fieldIvs[params_.timName.value()], timArr.data()),
        "nc_put_var_int (tim)");

  // Write the fields into the file
  // ------------------------------
  for (auto& field : fields) {
    // Get number of levels for this field
    const int nLevField = field.shape(1);

    // Create a rank 2 view of the field
    const auto fieldView = atlas::array::make_view<double, 2>(field);

    // Vector to hold the packed field
    std::vector<double> values(nLat*nLon*nLevField);

    // Loop over dimensions and pack the field
    for (size_t k = 0; k < nLevField; ++k) {
      for (size_t j = 0; j < nLat; ++j) {
        for (size_t i = 0; i < nLon; ++i) {
          values[k*nLat*nLon + j*nLon + i] = fieldView((nLat - 1 - j) * nLon + i, k);
        }
      }
    }

    // Write the field to the file
    nc_rc(nc_put_var_double(fileId, fieldIvs[field.name()], values.data()),
          "nc_put_var_double " + field.name());
  }

  // Close netCDF file
  // -----------------
  nc_rc(nc_close(fileId), "nc_close");
}

// -------------------------------------------------------------------------------------------------

void IOStructuredGrid::readStructuredFields(const std::string pathFile,
                                            atlas::FieldSet & fields,
                                            const util::DateTime & time,
                                            const eckit::LocalConfiguration & ioNames,
                                            const eckit::LocalConfiguration & ioScaling) const {
  // NetCDF IDs
  int fileId;

  // Open a file to read fields from
  // -------------------------------
  nc_rc(nc_open(pathFile.c_str(), NC_NOWRITE, &fileId), "nc_open " + pathFile);

  // Get file number of dimensions + their IDs
  // -----------------------------------------
  int ndims;
  nc_rc(nc_inq_ndims(fileId, &ndims), "nc_inq_ndims");

  std::vector<int> dimids(ndims);
  nc_rc(nc_inq_dimids(fileId, &ndims, dimids.data(), 0), "nc_inq_dimids");

  // Create regular grid for determining lat/lon values
  // --------------------------------------------------
  const atlas::RegularGrid regGrid(writeFunctionSpace_->grid());

  // Get grid dimensions
  // -------------------
  const int nLat = regGrid.ny();
  const int nLon = regGrid.nx();
  const int nLev = geom_.npz();
  const int nEdg = geom_.npz() + 1;
  const int nTim = 1;

  // Ensure that the lat and lon dimensions are found and have the correct lengths
  size_t dimSize;
  bool hasLat = false;
  bool hasLon = false;
  bool hasLev = false;
  bool hasEdg = false;
  bool hasTim = false;
  int latId;
  int lonId;
  int levId;
  int edgId;
  int timId;
  for (int i = 0; i < ndims; ++i) {
    // Get the name and size of the dimension
    char dimName[NC_MAX_NAME + 1];
    nc_rc(nc_inq_dim(fileId, dimids[i], dimName, &dimSize), "nc_inq_dim");

    if (std::string(dimName) == params_.latName.value().c_str()) {
      hasLat = true;
      latId = dimids[i];
      ASSERT(dimSize == nLat);
    } else if (std::string(dimName) == params_.lonName.value().c_str()) {
      hasLon = true;
      lonId = dimids[i];
      ASSERT(dimSize == nLon);
    } else if (std::string(dimName) == params_.levName.value().c_str()) {
      hasLev = true;
      levId = dimids[i];
      ASSERT(dimSize == nLev);
    } else if (std::string(dimName) == params_.edgName.value().c_str()) {
      hasEdg = true;
      edgId = dimids[i];
      ASSERT(dimSize == nEdg);
    } else if (std::string(dimName) == params_.timName.value().c_str()) {
      hasTim = true;
      timId = dimids[i];
      ASSERT(dimSize == nTim);
    }
  }

  // Ensure required dimensions were found
  ASSERT(hasLat);
  ASSERT(hasLon);
  ASSERT(hasLev);
  ASSERT(hasEdg);
  ASSERT(hasTim);

  // Read the fields from the file
  // -------------------------------
  for (auto & field : fields) {
    // Get IO name for this field
    std::string fieldName = field.name();
    if (ioNames.has(fieldName)) {
      fieldName = ioNames.getString(field.name());
    }

    // Get the variable ID for this field
    int varId;
    nc_rc(nc_inq_varid(fileId, fieldName.c_str(), &varId), "nc_inq_varid " + fieldName);

    // Get number of dimensions + their IDs
    int ndims;
    int dimids[NC_MAX_VAR_DIMS];
    nc_rc(nc_inq_var(fileId, varId,
                     nullptr,   // var name (unused)
                     nullptr,   // type (unused)
                     &ndims,
                     dimids,
                     nullptr),  // attributes (unused)
          "nc_inq_var");

    // Ensure that the field has either 3 or 4 dimensions
    ASSERT(ndims == 3 || ndims == 4);

    // Ensure that the dimensions are in the expected order
    size_t nLevField = 0;
    if ( ndims == 3 ) {
      ASSERT(dimids[0] == timId &&
             dimids[1] == latId &&
             dimids[2] == lonId);
      nLevField = 1;
    } else if ( ndims == 4 ) {
      ASSERT((dimids[0] == timId &&
              dimids[1] == levId &&
              dimids[2] == latId &&
              dimids[3] == lonId) ||
             (dimids[0] == timId &&
              dimids[1] == edgId &&
              dimids[2] == latId &&
              dimids[3] == lonId));

      nLevField = field.shape(1);
      if ( dimids[1] == edgId ) {
        ASSERT(nLevField == nEdg);
      } else {
        ASSERT(nLevField == nLev);
      }
    }

    // Read the variable data
    std::vector<double> values(field.size());

    nc_rc(nc_get_var_double(fileId, varId, values.data()), "nc_get_var_double " + fieldName);

    // Create field and unpack data into it
    auto fieldView = atlas::array::make_view<double, 2>(field);

    for (size_t k = 0; k < nLevField; ++k) {
      for (size_t j = 0; j < nLat; ++j) {
        for (size_t i = 0; i < nLon; ++i) {
          fieldView((nLat - 1 - j) * nLon + i, k) = values[ k*nLat*nLon + j*nLon + i ];
        }
      }
    }

    // Close file
    nc_rc(nc_close(fileId), "nc_close");
  }
}

// -------------------------------------------------------------------------------------------------

void IOStructuredGrid::print(std::ostream & os) const {
  os << classname() << " IO using Atlas Structured Grid";
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
