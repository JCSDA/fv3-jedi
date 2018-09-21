/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <mpi.h>
#include <unistd.h>
#include <string>
#include "eckit/config/Configuration.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "src/Utilities/UtilitiesFV3JEDI.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {

// -----------------------------------------------------------------------------

void stageFv3Files(const eckit::Configuration &conf) {
  oops::Log::trace() << "Staging fv3 input.nml and field_table" << std::endl;

  // Get processor ID
  int world_rank = oops::mpi::comm().rank();

  // Only one processor needs to move the files
  if (world_rank == 0) {
    // Remove anything currently present
    std::remove("input.nml");
    std::remove("field_table");

    // User provided input files for this geom/state/model etc
    if (conf.has("nml_file")) {
      std::string nml_file = conf.getString("nml_file");
      symlink(nml_file.c_str(), "./input.nml");
    } else {
      ABORT("input.nml not in configuration");
    }
    if (conf.has("trc_file")) {
      std::string trc_file = conf.getString("trc_file");
      symlink(trc_file.c_str(), "./field_table");
    } else {
      ABORT("field_table not in configuration");
    }

    // User may also be requesting the tlm/adm nml file
    if (conf.has("nml_file_pert")) {
      oops::Log::trace() << "Also staging fv3 inputpert.nml" << std::endl;
      std::remove("inputpert.nml");
      std::string nml_file_pert = conf.getString("nml_file_pert");
      symlink(nml_file_pert.c_str(), "./inputpert.nml");
    }
  }

  // Nobody moves until files are in place
  oops::mpi::comm().barrier();
}

// -----------------------------------------------------------------------------

void stageFv3Input(const eckit::Configuration &conf) {
  oops::Log::trace() << "Staging fv3 input.nml" << std::endl;

  // Get processor ID
  int world_rank = oops::mpi::comm().rank();

  // Only one processor needs to move the files
  if (world_rank == 0) {
    // Remove anything currently present
    std::remove("input.nml");

    // User provided input files for this geom/state/model etc
    if (conf.has("nml_file")) {
      std::string nml_file = conf.getString("nml_file");
      symlink(nml_file.c_str(), "./input.nml");
    } else {
      ABORT("input.nml not in configuration");
    }
  }

  // Nobody moves until files are in place
  oops::mpi::comm().barrier();
}

// -----------------------------------------------------------------------------

void removeFv3Files() {
  oops::Log::trace() << "Removing fv3 input.nml and field_table" << std::endl;

  // Get processor ID
  int world_rank = oops::mpi::comm().rank();

  // No file deletion until everyone catches up
  oops::mpi::comm().barrier();

  // Only one processor needs to move the files
  if (world_rank == 0) {
    std::remove("input.nml");
    std::remove("field_table");

    // And inputpert if it is there
    std::remove("inputpert.nml");
  }
}

// -----------------------------------------------------------------------------

void removeFv3Input() {
  oops::Log::trace() << "Removing fv3 input.nml" << std::endl;

  // Get processor ID
  int world_rank = oops::mpi::comm().rank();

  // No file deletion until everyone catches up
  oops::mpi::comm().barrier();

  // Only one processor needs to move the files
  if (world_rank == 0) {
    std::remove("input.nml");
  }
}

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
