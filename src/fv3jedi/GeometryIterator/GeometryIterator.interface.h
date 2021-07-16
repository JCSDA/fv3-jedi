/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_GEOMETRYITERATOR_GEOMETRYITERATOR_INTERFACE_H_
#define FV3JEDI_GEOMETRYITERATOR_GEOMETRYITERATOR_INTERFACE_H_

#include "fv3jedi/Utilities/interface.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace fv3jedi {

extern "C" {

  void fv3jedi_geom_iter_setup_f90(F90iter &, const F90geom &,
                                const int &, const int &);
  void fv3jedi_geom_iter_clone_f90(F90iter &, const F90iter &);
  void fv3jedi_geom_iter_delete_f90(F90iter &);
  void fv3jedi_geom_iter_equals_f90(const F90iter &, const F90iter&, int &);
  void fv3jedi_geom_iter_current_f90(const F90iter &, double &, double &);
  void fv3jedi_geom_iter_next_f90(const F90iter &);
  void fv3jedi_geom_iter_orography_f90(const F90iter &, double &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_GEOMETRYITERATOR_GEOMETRYITERATOR_INTERFACE_H_
