
#include "RunFV3JEDI.h"

#include "Fortran.h"
#include "util/Logger.h"
#include "oops/runs/Run.h"
#include "eckit/config/Configuration.h"
#include <fstream>

namespace fv3jedi {

// -----------------------------------------------------------------------------

RunFV3JEDI::RunFV3JEDI(int argc, char ** argv) : oops::Run(argc, argv) {
  oops::Log::trace() << "Creating RunFV3JEDI" << std::endl;
  const eckit::Configuration * conf = &config();

  std::remove("input.nml");
  std::string name = "input.nml";
  std::ofstream f2(name);

  fv3jedi_setup_f(&conf);

  std::remove("input.nml");

  oops::Log::trace() << "RunFV3JEDI created" << std::endl;
}

// -----------------------------------------------------------------------------

RunFV3JEDI::~RunFV3JEDI() {
  oops::Log::trace() << "Destructing RunFV3JEDI" << std::endl;
 
  std::string name = "input.nml";
  std::ofstream f2(name);

  fv3jedi_finalize_f();

  std::remove("input.nml");

  oops::Log::trace() << "MPI finalized, RunFV3JEDI destructed" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
