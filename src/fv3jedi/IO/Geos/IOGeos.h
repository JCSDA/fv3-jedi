/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "oops/util/DateTime.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "fv3jedi/IO/Utils/IOBase.h"
#include "IOGeos.interface.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

class IOGeosParameters : public IOParametersBase {
  OOPS_CONCRETE_PARAMETERS(IOGeosParameters, IOParametersBase)

 public:
  // Data path for files being read
  oops::Parameter<std::string> datapath{"datapath", "path to location of files to be read",
                                        "./", this};

  // Option to clobber existing files
  oops::Parameter<bool> clobber{"clobber", "clobber existing netcdf files", true, this};

  // Whether the tile is a dimension in the file
  oops::Parameter<bool> tiledim{"tiledim", "file uses the tile dimension", true, this};

  // Filenames to be read
  oops::Parameter<std::string> filename_bkgd{"filename_bkgd", "filename_bkgd",
                                             "bkg", this};
  oops::Parameter<std::string> filename_crtm{"filename_crtm", "filename_crtm",
                                             "crtmsrf", this};
  oops::Parameter<std::string> filename_core{"filename_core", "filename_core",
                                             "fvcore_internal_rst", this};
  oops::Parameter<std::string> filename_mois{"filename_mois", "filename_mois",
                                             "moist_internal_rst", this};
  oops::Parameter<std::string> filename_surf{"filename_surf", "filename_surf",
                                             "surf_import_rst", this};

  // Input filename may be templated with datetimes
  oops::Parameter<bool> filename_is_datetime_templated{"filename is datetime templated",
                                                       "filename is datetime templated",
                                                       false, this};

  // Force tell the system that surface pressure is in the file
  oops::Parameter<bool> psinfile{"psinfile",
                                 "tell the system surface pressure is in the file",
                                 false, this};

  // Optionally the config may contain member
  oops::OptionalParameter<int> member{"member", "ensemble member number", this};

  // Things BUMP puts in write config not needed in IO but specified to avoid failures
  oops::OptionalParameter<std::string> bumpparameter{"parameter", "bump parameter", this};
  oops::OptionalParameter<util::DateTime> bumpdatetime{"date", "bump datetime", this};
};

// -------------------------------------------------------------------------------------------------
class IOGeos : public IOBase, private util::ObjectCounter<IOGeos> {
 public:
  static const std::string classname() {return "fv3jedi::IOGeos";}

  typedef IOGeosParameters Parameters_;

  IOGeos(const Geometry &, const Parameters_ &);
  ~IOGeos();
  void read(State &) const override;
  void read(Increment &) const override;
  void write(const State &) const override;
  void write(const Increment &) const override;

 private:
  F90iogeos objectKeyForFortran_;
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
