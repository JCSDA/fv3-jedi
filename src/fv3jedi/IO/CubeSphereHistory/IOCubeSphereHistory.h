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

  // Optionally the config may contain member
  oops::OptionalParameter<int> member{"member", "ensemble member number", this};

  // Things BUMP puts in write config not needed in IO but specified to avoid failures
  oops::OptionalParameter<std::string> bumpparameter{"parameter", "bump parameter", this};
  oops::OptionalParameter<util::DateTime> bumpdatetime{"date", "bump datetime", this};
};

// -------------------------------------------------------------------------------------------------
class IOCubeSphereHistory : public IOBase, private util::ObjectCounter<IOCubeSphereHistory> {
 public:
  static const std::string classname() {return "fv3jedi::IOCubeSphereHistory";}

  typedef IOCubeSphereHistoryParameters Parameters_;

  IOCubeSphereHistory(const Geometry &, const Parameters_ &);
  ~IOCubeSphereHistory();
  void read(State &) const override;
  void read(Increment &) const override;
  void write(const State &) const override;
  void write(const Increment &) const override;

 private:
  F90IOCubeSphereHistory objectKeyForFortran_;
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
