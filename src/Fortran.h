/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3_JEDI_SRC_FORTRAN_H_
#define FV3_JEDI_SRC_FORTRAN_H_

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
  class Duration;
}

namespace fv3jedi {

// Geometry key type
typedef int F90geom;
// Model key type
typedef int F90model;
// Locations key type
typedef int F90locs;
// Goms key type
typedef int F90goms;
// Fields key type
typedef int F90flds;
// Trajectory key type
typedef int F90traj;
// Background error covariance key type
typedef int F90bmat;
// Localization matrix
typedef int F90lclz;
// ObOp trajectory
typedef int F90ootrj;
// VarChange key
typedef int F90vcha;


/// Interface to Fortran FV3JEDI model
/*!
 * The core of the FV3JEDI model is coded in Fortran.
 * Here we define the interfaces to the Fortran code.
 */

extern "C" {
// -----------------------------------------------------------------------------
//  For all FV3 runs
// -----------------------------------------------------------------------------
  void fv3jedi_setup_f(const eckit::Configuration * const *);
  void fv3jedi_finalize_f();

// -----------------------------------------------------------------------------
//  Geometry
// -----------------------------------------------------------------------------
  void fv3jedi_geo_setup_f90(F90geom &, const eckit::Configuration * const *);
  void fv3jedi_geo_clone_f90(const F90geom &, F90geom &);
  void fv3jedi_geo_info_f90(const F90geom &);
  // void fv3jedi_geo_info_f90(const F90geom &, int &, int &);
  void fv3jedi_geo_delete_f90(F90geom &);

// -----------------------------------------------------------------------------
//  Model
// -----------------------------------------------------------------------------
  void fv3jedi_model_setup_f90(const eckit::Configuration * const *,
                                const F90geom &, F90model &);
  void fv3jedi_model_delete_f90(F90model &);

  void fv3jedi_model_prepare_integration_f90(const F90model &, const F90flds &);
  void fv3jedi_model_prepare_integration_tl_f90(const F90model &,
                                                 const F90flds &);
  void fv3jedi_model_prepare_integration_ad_f90(const F90model &,
                                                 const F90flds &);

  void fv3jedi_model_propagate_f90(const F90model &, const F90flds &);
  void fv3jedi_model_prop_traj_f90(const F90model &, const F90flds &,
                                    F90traj &);
  void fv3jedi_model_propagate_tl_f90(const F90model &, const F90flds &,
                                       const F90traj &);
  void fv3jedi_model_propagate_ad_f90(const F90model &, const F90flds &,
                                       const F90traj &);

  void fv3jedi_model_wipe_traj_f90(F90traj &);
  void fv3jedi_traj_minmaxrms_f90(const F90traj &, double &);

// -----------------------------------------------------------------------------
//  Fields
// -----------------------------------------------------------------------------
  void fv3jedi_field_create_f90(F90flds &, const F90geom &,
                                  const eckit::Configuration * const *);
  void fv3jedi_field_delete_f90(F90flds &);

  void fv3jedi_field_copy_f90(const F90flds &, const F90flds &);
  void fv3jedi_field_zero_f90(const F90flds &);
  void fv3jedi_field_self_add_f90(const F90flds &, const F90flds &);
  void fv3jedi_field_self_sub_f90(const F90flds &, const F90flds &);
  void fv3jedi_field_self_mul_f90(const F90flds &, const double &);
  void fv3jedi_field_axpy_f90(const F90flds &, const double &, const F90flds &);
  void fv3jedi_field_dot_prod_f90(const F90flds &, const F90flds &, double &);
  void fv3jedi_field_self_schur_f90(const F90flds &, const F90flds &);
  void fv3jedi_field_random_f90(const F90flds &);

  void fv3jedi_field_add_incr_f90(const F90flds &, const F90flds &);
  void fv3jedi_field_diff_incr_f90(const F90flds &, const F90flds &,
                                    const F90flds &);

  void fv3jedi_field_change_resol_f90(const F90flds &, const F90flds &);

  void fv3jedi_field_read_file_f90(const F90flds &,
                                    const eckit::Configuration * const *,
                                    util::DateTime * const *);
  void fv3jedi_field_analytic_init_f90(const F90flds &, const F90geom &,
                                        const eckit::Configuration * const *,
                                        util::DateTime * const *);
  void fv3jedi_field_write_file_f90(const F90flds &,
                                     const eckit::Configuration * const *,
                                     const util::DateTime * const *);

  void fv3jedi_field_getvalues_notraj_f90(const F90flds &, const F90locs &,
                        const eckit::Configuration * const *, const F90goms &);
  void fv3jedi_field_getvalues_f90(const F90flds &, const F90locs &,
                        const eckit::Configuration * const *, const F90goms &,
                        const F90ootrj &);
  void fv3jedi_field_getvalues_tl_f90(const F90flds &, const F90locs &,
                        const eckit::Configuration * const *, const F90goms &,
                        const F90ootrj &);
  void fv3jedi_field_getvalues_ad_f90(const F90flds &, const F90locs &,
                        const eckit::Configuration * const *, const F90goms &,
                        const F90ootrj &);

  void fv3jedi_field_gpnorm_f90(const F90flds &, const int &, double &);
  void fv3jedi_field_sizes_f90(const F90flds &, int &, int &, int &);
  void fv3jedi_field_rms_f90(const F90flds &, double &);
  void fv3jedi_field_ug_coord_f90(const F90flds &, const int &);
  void fv3jedi_field_field_to_ug_f90(const F90flds &, const int &);
  void fv3jedi_field_field_from_ug_f90(const F90flds &, const int &);

  void fv3jedi_field_dirac_f90(const F90flds &,
                                const eckit::Configuration * const *);

  void fv3jedi_getvaltraj_setup_f90(const F90ootrj &);
  void fv3jedi_getvaltraj_delete_f90(const F90ootrj &);

// -----------------------------------------------------------------------------
//  Background error
// -----------------------------------------------------------------------------
  void fv3jedi_b_setup_f90(F90bmat &, const eckit::Configuration * const *,
                             const F90geom &);
  void fv3jedi_b_delete_f90(F90bmat &);

  void fv3jedi_b_linearize_f90(const F90bmat &,
                                 const eckit::Configuration * const *);

  void fv3jedi_b_mult_f90(const F90bmat &, const F90flds &, const F90flds &);
  void fv3jedi_b_invmult_f90(const F90bmat &, const F90flds &, const F90flds &);

  void fv3jedi_b_randomize_f90(const F90bmat &, const F90flds &);

// -----------------------------------------------------------------------------
//  Change variable for B matrix
// -----------------------------------------------------------------------------

  void fv3jedi_varchange_setup_f90(const F90vcha &,
                                   const eckit::Configuration * const *);
  void fv3jedi_varchange_delete_f90(F90vcha &);
  void fv3jedi_varchange_linearize_f90(F90vcha &,const F90geom &,
                                       const F90flds &);
  void fv3jedi_varchange_transform_f90(const F90vcha &, const F90flds &,
                                       const F90flds &);
  void fv3jedi_varchange_transformadjoint_f90(const F90vcha &, const F90flds &,
                                       const F90flds &);
  void fv3jedi_varchange_transforminverse_f90(const F90vcha &, const F90flds &,
                                       const F90flds &);
  void fv3jedi_varchange_transforminverseadjoint_f90(const F90vcha &,
                                              const F90flds &, const F90flds &);

// -----------------------------------------------------------------------------
//  Localization matrix
// -----------------------------------------------------------------------------
  void fv3jedi_localization_setup_f90(F90lclz &,
                                 const eckit::Configuration * const *,
                                 const F90geom &);
  void fv3jedi_localization_delete_f90(F90lclz &);
  void fv3jedi_localization_mult_f90(const F90lclz &, const F90flds &);

}  // extern "C"
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3_JEDI_SRC_FORTRAN_H_
