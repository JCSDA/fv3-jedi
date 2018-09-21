/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <fstream>

#include "eckit/config/Configuration.h"
#include "oops/runs/Run.h"
#include "oops/util/Logger.h"

#include "src/Run/RunFV3JEDI.h"
#include "RunFV3JEDIFortran.h"
#include "UtilitiesFV3JEDI.h"

namespace fv3jedi {

// -----------------------------------------------------------------------------

RunFV3JEDI::RunFV3JEDI(int argc, char ** argv) : oops::Run(argc, argv) {
  oops::Log::trace() << "Creating RunFV3JEDI" << std::endl;
  const eckit::Configuration * conf = &config();

  stageFv3Input(config());
  fv3jedi_setup_f(&conf);
  removeFv3Input();

  oops::Log::trace() << "RunFV3JEDI created" << std::endl;
}

// -----------------------------------------------------------------------------

RunFV3JEDI::~RunFV3JEDI() {
  oops::Log::trace() << "Destructing RunFV3JEDI" << std::endl;
  fv3jedi_finalize_f();
  oops::Log::trace() << "MPI finalized, RunFV3JEDI destructed" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
