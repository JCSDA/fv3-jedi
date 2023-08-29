/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

#include "fv3jedi/FieldMetadata/FieldsMetadata.h"
#include "fv3jedi/Utilities/interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace fv3jedi {

extern "C" {

  void fv3jedi_geom_initialize_f90(const eckit::LocalConfiguration &, const eckit::mpi::Comm *);
  void fv3jedi_geom_setup_f90(F90geom &, const eckit::Configuration &,
                             const eckit::mpi::Comm *, int &);
  void fv3jedi_geom_addfmd_f90(F90geom &, FieldsMetadata *);
  void fv3jedi_geom_set_lonlat_f90(const F90geom &, atlas::field::FieldSetImpl *,
                                   const bool &);
  void fv3jedi_geom_set_functionspace_pointer_f90(const F90geom &,
                                                  atlas::functionspace::FunctionSpaceImpl *,
                                                  atlas::functionspace::FunctionSpaceImpl *);
  void fv3jedi_geom_set_and_fill_geometry_fields_f90(const F90geom &, atlas::field::FieldSetImpl *);
  void fv3jedi_geom_clone_f90(F90geom &, const F90geom &, const FieldsMetadata *);
  void fv3jedi_geom_print_f90(const F90geom &, int &);
  void fv3jedi_geom_delete_f90(F90geom &);
  void fv3jedi_geom_start_end_f90(const F90geom &, int &, int &, int &, int &, int &,
                                  int &, int &);
  void fv3jedi_geom_verticalCoord_f90(const F90geom &, double &, int &, double &);
  int fv3jedi_geom_iterator_dimension_f90(const F90geom &, int &);
  void fv3jedi_geom_get_data_f90(const F90geom &, const int &, double *, double *, double &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
