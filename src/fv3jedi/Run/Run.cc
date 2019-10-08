/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <fstream>
#include <string>

#include "eckit/config/Configuration.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Run.h"
#include "oops/util/Logger.h"

#include "fv3jedi/Run/Run.h"
#include "fv3jedi/Utilities/Utilities.h"

namespace fv3jedi {

// -----------------------------------------------------------------------------

Run::Run(int argc, char ** argv) : oops::Run(argc, argv) {
  oops::Log::trace() << "Creating Run" << std::endl;
  const eckit::Configuration * conf = &config();

  /* NOTE: the following lines should move to another place
  after the EDA has been constructed and the communicators exist.
  Right now we use MPI_COMM_WORLD */

  stageFv3Files(config());
  std::string commName = oops::mpi::comm().name();
  fv3jedi_setup_f(&conf, commName.size(), commName.c_str());
  removeFv3Files();

  oops::Log::trace() << "Run created" << std::endl;
}

// -----------------------------------------------------------------------------

Run::~Run() {
  oops::Log::trace() << "Destructing Run" << std::endl;
  fv3jedi_finalize_f();
  oops::Log::trace() << "MPI finalized, Run destructed" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
