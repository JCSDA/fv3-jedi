/*
 * (C) Copyright 2021-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace fv3jedi {
  typedef int F90IOCubeSphereHistory;
  extern "C" {
  void fv3jedi_io_cube_sphere_history_create_f90(F90IOCubeSphereHistory &,
                                                 const eckit::Configuration &,
                                                 const F90geom &);
  void fv3jedi_io_cube_sphere_history_delete_f90(F90IOCubeSphereHistory &);
  void fv3jedi_io_cube_sphere_history_read_state_f90(const F90IOCubeSphereHistory &, F90state &);
  void fv3jedi_io_cube_sphere_history_read_increment_f90(const F90IOCubeSphereHistory &, F90inc &);
  void fv3jedi_io_cube_sphere_history_write_state_f90(const F90IOCubeSphereHistory &,
                                                      const F90state &);
  void fv3jedi_io_cube_sphere_history_write_increment_f90(const F90IOCubeSphereHistory &,
                                                          const F90inc &);
  }  // extern "C"
}  // namespace fv3jedi
