/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "oops/util/DateTime.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "fv3jedi/IO/Utils/IOBase.h"
#include "IOCubeSphereHistory.interface.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

class IOCubeSphereHistoryParameters : public IOParametersBase {
  OOPS_CONCRETE_PARAMETERS(IOCubeSphereHistoryParameters, IOParametersBase)

 public:
  // Names of files to be read/written to
  oops::Parameter<std::string> provider{"provider", "name of the model", "geos", this};

  // Filenames provided as a list
  oops::OptionalParameter<std::vector<std::string>> filenames{"filenames",
                                                              "names of the files to be read",
                                                              this};

  // Single filename provided
  oops::OptionalParameter<std::string> filename{"filename",
                                                "name of the file to be read",
                                                this};

  // Path prepended to all files
  oops::Parameter<std::string> datapath{"datapath", "path to location of files to be read",
                                        "./", this};

  // Option to clobber existing files
  oops::OptionalParameter<std::vector<bool>> clobber{"clobber existing files",
                                                     "clobber existing files", this};

  // Whether the tile is a dimension in the file
  oops::OptionalParameter<std::vector<bool>> tiledim{"tile is a dimension",
                                             "tile is a dimension", this};

  // Name of the X Dimension in the file
  oops::OptionalParameter<std::vector<std::string>> xdim{"x dimension name",
                                                         "x dimension name",
                                                         this};

  // Name of the Y Dimension in the file
  oops::OptionalParameter<std::vector<std::string>> ydim{"y dimension name",
                                                         "y dimension name",
                                                         this};

  // Name of the Z Full Dimension in the file
  oops::OptionalParameter<std::vector<std::string>> zfdim{"z full dimension name",
                                                          "z full dimension name",
                                                          this};

  // Name of the Z Half Dimension in the file
  oops::OptionalParameter<std::vector<std::string>> zhdim{"z half dimension name",
                                                          "z half dimension name",
                                                          this};

  // Set date/time on read
  oops::OptionalParameter<bool> setDateTime{"set datetime on read", "set datetime on read", this};

  // Optionally the config may contain member
  oops::OptionalParameter<int> member{"member", "ensemble member number", this};

  // Optional list of fields to write out
  oops::OptionalParameter<std::vector<std::string>> fieldsToWrite{"fields to write",
                                                                  "names of the fields to write",
                                                                  this};

  // Floating point precision in bytes for NetCDF write
  oops::OptionalParameter<int> floatPrecision{"float precision in bytes",
                                              "number of bytes of floating point precision",
                                              this};

  // Compute pressure at the edges from pressure at the surface (instead of reading it)
  oops::OptionalParameter<bool> computeP{"compute edge pressure from surface pressure",
                                         "compute edge pressure from surface pressure",
                                         this};

  // Maximum allowable difference in the Geometry lat/lon compared to the file lat/lon
  // In practice users should expect differences order 1e-12 or smaller if everything
  // is in double precision. In practice models may produce files at lower precision.
  // Differences smaller than 1e-6 should be sufficient to assess that the geometry of
  // the model producing the file being read is the same at the one in fv3-jedi.
  oops::Parameter<double> maxDiff{"max allowable geometry difference",
                                  "max allowable geometry difference", 1e-6, this};
};

// -------------------------------------------------------------------------------------------------
class IOCubeSphereHistory : public IOBase, private util::ObjectCounter<IOCubeSphereHistory> {
 public:
  static const std::string classname() {return "fv3jedi::IOCubeSphereHistory";}

  typedef IOCubeSphereHistoryParameters Parameters_;

  IOCubeSphereHistory(const Geometry &, const Parameters_ &);
  ~IOCubeSphereHistory();
  void read(State &, const eckit::LocalConfiguration &,
            const eckit::LocalConfiguration &) const override;
  void read(Increment &, const eckit::LocalConfiguration &,
            const eckit::LocalConfiguration &) const override;
  void write(const State &, const eckit::LocalConfiguration &,
             const eckit::LocalConfiguration &) const override;
  void write(const Increment &, const eckit::LocalConfiguration &,
             const eckit::LocalConfiguration &) const override;

 private:
  F90IOCubeSphereHistory objectKeyForFortran_;
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
