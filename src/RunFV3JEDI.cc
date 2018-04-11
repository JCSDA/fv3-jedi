
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

  // Setup an input.nml
  // This should be replaced with symbolic link to file 
  std::remove("input.nml");
  std::ofstream nmlfile("input.nml");
  std::string nmlstr;
  nmlstr = "&fms_io_nml\n"                             ; nmlfile << nmlstr;
  nmlstr = "       checksum_required       = F\n"      ; nmlfile << nmlstr;
  nmlstr = "/\n"                                       ; nmlfile << nmlstr;
  nmlstr = "\n"                                        ; nmlfile << nmlstr;
  nmlstr = "&fms_nml\n"                                ; nmlfile << nmlstr;
  nmlstr = "       print_memory_usage=.false.\n"       ; nmlfile << nmlstr;
  nmlstr = "       domains_stack_size = 24000000\n"    ; nmlfile << nmlstr;
  nmlstr = "/\n"                                       ; nmlfile << nmlstr;
  nmlstr = "\n"                                        ; nmlfile << nmlstr;  
  nmlstr = "&mpp_io_nml \n"                            ; nmlfile << nmlstr;
  nmlstr = "       header_buffer_val       = 16384 \n" ; nmlfile << nmlstr;
  nmlstr = "       global_field_on_root_pe = T \n"     ; nmlfile << nmlstr;
  nmlstr = "       io_clocks_on            = F \n"     ; nmlfile << nmlstr;
  nmlstr = "       shuffle                 = 0 \n"     ; nmlfile << nmlstr;
  nmlstr = "       deflate_level           = -1 \n"    ; nmlfile << nmlstr;
  nmlstr = "       cf_compliance           = F \n"     ; nmlfile << nmlstr;
  nmlstr = "/\n"; nmlfile << nmlstr;
  nmlfile.close();

  fv3jedi_setup_f(&conf);

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
