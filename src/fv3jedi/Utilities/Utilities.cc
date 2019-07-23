/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <mpi.h>
#include <unistd.h>

#include <fstream>
#include <iostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/parallel/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "fv3jedi/Utilities/Utilities.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {

// -----------------------------------------------------------------------------

void stageFv3Files(const eckit::Configuration &conf) {
  oops::Log::trace() << "Staging files for FV3" << std::endl;

  // Get processor ID
  int world_rank = oops::mpi::comm().rank();

  // Only one processor needs to move the files
  if (world_rank == 0) {
    // User provided input files for this geom/state/model etc
    delete_file("input.nml");
    if (conf.has("nml_file")) {
      oops::Log::debug() << "Staging input.nml" << std::endl;
      std::string nml_file = conf.getString("nml_file");
      symlink(nml_file.c_str(), "./input.nml");
    } else {
      ABORT("nml_file not in configuration");
    }

    // User may also be requesting the field_table to be staged
    delete_file("field_table");
    if (conf.has("trc_file")) {
      oops::Log::debug() << "Staging field_table" << std::endl;
      std::string trc_file = conf.getString("trc_file");
      symlink(trc_file.c_str(), "./field_table");
    }

    // User may also be requesting the tlm/adm nml file
    delete_file("inputpert.nml");
    if (conf.has("nml_file_pert")) {
      oops::Log::debug() << "Staging inputpert.nml" << std::endl;
      std::string nml_file_pert = conf.getString("nml_file_pert");
      symlink(nml_file_pert.c_str(), "./inputpert.nml");
    }
  }

  // Nobody moves until files are in place
  oops::mpi::comm().barrier();
}

// -----------------------------------------------------------------------------

void removeFv3Files() {
  oops::Log::trace() << "Removing staged fv3 files" << std::endl;

  // Get processor ID
  int world_rank = oops::mpi::comm().rank();

  // No file deletion until everyone catches up
  oops::mpi::comm().barrier();

  // Only one processor needs to move the files
  if (world_rank == 0) {
    delete_file("input.nml");
    delete_file("field_table");
    delete_file("inputpert.nml");
  }
}

// -----------------------------------------------------------------------------

void delete_file(const char *fileName)
{
  std::ifstream infile(fileName);
  if (infile.good()) {
    oops::Log::debug() << "Removing: " << fileName << std::endl;
    std::remove(fileName);
  }
}

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
