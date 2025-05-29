/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <netcdf.h>

#include <map>
#include <vector>

#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "oops/base/GeometryData.h"
#include "oops/util/Logger.h"
#include "oops/util/stringFunctions.h"
#include "oops/util/Timer.h"

#include "fv3jedi/FieldMetadata/FieldsMetadata.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/IO/StructuredGrid/IOStructuredGrid.h"
#include "fv3jedi/State/State.h"

namespace fv3jedi {
// -------------------------------------------------------------------------------------------------
static IOMaker<IOStructuredGrid> makerIOStructuredGrid_("structured grid");
static IOMaker<IOStructuredGrid> makerIOAuxGrid_("auxgrid");
// -------------------------------------------------------------------------------------------------
IOStructuredGrid::IOStructuredGrid(const Geometry & geom, const Parameters_ & params)
  : IOBase(geom, params.toConfiguration()), interpolator_(), params_(params), gridStr_(""),
    geom_(geom), writeFunctionSpace_() {
  util::Timer timer(classname(), "IOStructuredGrid");
  oops::Log::trace() << classname() << " constructor starting" << std::endl;

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

  // Create a generic interpolator for converting to the structured grid
  // -------------------------------------------------------------------
  interpolator_.reset(new oops::GlobalInterpolator(params.toConfiguration(), geomData,
                                                   *writeFunctionSpace_,
                                                   geom.getComm()));
  oops::Log::trace() << classname() << " constructor done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
IOStructuredGrid::~IOStructuredGrid() {
  util::Timer timer(classname(), "~IOStructuredGrid");
  oops::Log::trace() << classname() << " destructor starting" << std::endl;
  oops::Log::trace() << classname() << " destructor done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void IOStructuredGrid::read(State & x, const eckit::LocalConfiguration & fileionames,
                                const eckit::LocalConfiguration & fileioscaling) const {
  ABORT("IOStructuredGrid::read(State) not implemented");
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

  // Write to disk if rank 0
  if (geom_.getComm().rank() == 0) {
    this->writeStructuredFields(fieldsGeographic, obj.validTime(), fileionames, fileioscaling);
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

static inline void nc_rc(const int return_code, const std::string & operation) {
  // If there was a failure of the netCDF operation, abort with the error message
  if (return_code) {
    ABORT("IOStructuredGrid netCDF operation \'" + operation + "\' failed with error: "
          + nc_strerror(return_code));
  }
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
    latArr[i] = regGrid.y(i);
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

void IOStructuredGrid::print(std::ostream & os) const {
  os << classname() << " IO using Atlas Structured Grid";
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
