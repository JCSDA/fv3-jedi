/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_GEOMETRY_GEOMETRY_INTERFACE_H_
#define FV3JEDI_GEOMETRY_GEOMETRY_INTERFACE_H_

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "eckit/mpi/Comm.h"
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

  void fv3jedi_geo_setup_f90(F90geom &, const eckit::Configuration * const *,
                             const eckit::mpi::Comm *);
  void fv3jedi_geo_create_atlas_grid_conf_f90(const F90geom &,
                                              const eckit::Configuration * const *);
  void fv3jedi_geo_set_atlas_functionspace_pointer_f90(const F90geom &,
                                                       atlas::functionspace::FunctionSpaceImpl *);
  void fv3jedi_geo_fill_atlas_fieldset_f90(const F90geom &,
                                           atlas::field::FieldSetImpl *);
  void fv3jedi_geo_set_atlas_fieldset_pointer_f90(const F90geom &,
                                                  atlas::field::FieldSetImpl *);
  void fv3jedi_geo_clone_f90(const F90geom &,
                             F90geom &);
  void fv3jedi_geo_info_f90(const F90geom &);
  void fv3jedi_geo_delete_f90(F90geom &);
  void fv3jedi_geo_start_end_f90(const F90geom &, int &, int &, int &, int &,
                                 int &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_GEOMETRY_GEOMETRY_INTERFACE_H_
