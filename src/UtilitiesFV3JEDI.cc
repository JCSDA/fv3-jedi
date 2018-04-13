/*
 * (C) Copyright 2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <mpi.h>
#include "eckit/config/Configuration.h"
#include "util/Logger.h"
#include "UtilitiesFV3JEDI.h"
#include <unistd.h>

namespace eckit {
  class Configuration;
}

namespace fv3jedi {

 void stageFv3Files (const eckit::Configuration &conf) {

   oops::Log::trace() << "Staging fv3 input.nml and field_table" << std::endl;

   // Get processor ID
   int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

   // Only one processor needs to move the files
   if (world_rank == 0) {

     // Remove anything currently present
     std::remove("input.nml");
     std::remove("field_table");

     // User provided input files for this geom/state/model etc
     std::string nml_file = conf.getString("nml_file");
     std::string trc_file = conf.getString("trc_file");

     // Create links in the cwd 
     symlink(nml_file.c_str(), "./input.nml");
     symlink(trc_file.c_str(), "./field_table");

   }

   // Nobody moves until files are in place
   MPI_Barrier(MPI_COMM_WORLD);

 }

 void removeFv3Files () {

   oops::Log::trace() << "Removing fv3 input.nml and field_table" << std::endl;

   // Get processor ID
   int world_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

   // Only one processor needs to move the files
   if (world_rank == 0) {
     std::remove("input.nml");
     std::remove("field_table");
   }

 }

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
