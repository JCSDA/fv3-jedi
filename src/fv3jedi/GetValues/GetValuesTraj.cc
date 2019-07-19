/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <mpi.h>

#include "eckit/config/Configuration.h"
#include "oops/util/Logger.h"

#include "fv3jedi/GetValues/GetValuesTraj.h"

// -----------------------------------------------------------------------------
namespace fv3jedi {
// -----------------------------------------------------------------------------
GetValuesTrajMatrix::GetValuesTrajMatrix() {
  oops::Log::trace() << "GetValuesTrajMatrix constructor starting"
                     << std::endl;
  fv3jedi_getvalues_traj_setup_f90(keyGetValuesTraj_);
  oops::Log::trace() << "GetValuesTrajMatrix constructor done"
                     << keyGetValuesTraj_ << std::endl;
}
// -----------------------------------------------------------------------------
GetValuesTrajMatrix::~GetValuesTrajMatrix() {
  oops::Log::trace() << "GetValuesTrajMatrix destructor starting"
                     << std::endl;
  fv3jedi_getvalues_traj_delete_f90(keyGetValuesTraj_);
  oops::Log::trace() << "GetValuesTrajMatrix destructor done" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
