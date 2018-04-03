/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "FieldsFV3JEDI.h"

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/generic/UnstructuredGrid.h"
#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "util/Logger.h"
#include "Fortran.h"
#include "GeometryFV3JEDI.h"
#include "util/DateTime.h"

// -----------------------------------------------------------------------------
namespace fv3jedi {
// -----------------------------------------------------------------------------
FieldsFV3JEDI::FieldsFV3JEDI(const GeometryFV3JEDI & geom, const oops::Variables & vars,
                           const util::DateTime & time):
  geom_(new GeometryFV3JEDI(geom)), vars_(vars), time_(time)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  fv3jedi_field_create_f90(keyFlds_, geom_->toFortran(), &conf);
}
// -----------------------------------------------------------------------------
FieldsFV3JEDI::FieldsFV3JEDI(const FieldsFV3JEDI & other, const bool copy)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  fv3jedi_field_create_f90(keyFlds_, geom_->toFortran(), &conf);
  if (copy) {
    fv3jedi_field_copy_f90(keyFlds_, other.keyFlds_);
  } else {
    fv3jedi_field_zero_f90(keyFlds_);
  }
}
// -----------------------------------------------------------------------------
FieldsFV3JEDI::FieldsFV3JEDI(const FieldsFV3JEDI & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  fv3jedi_field_create_f90(keyFlds_, geom_->toFortran(), &conf);
  fv3jedi_field_copy_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsFV3JEDI::FieldsFV3JEDI(const FieldsFV3JEDI & other, const GeometryFV3JEDI & geom)
  : geom_(new GeometryFV3JEDI(geom)), vars_(other.vars_), time_(other.time_)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  fv3jedi_field_create_f90(keyFlds_, geom_->toFortran(), &conf);
  fv3jedi_field_change_resol_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsFV3JEDI::FieldsFV3JEDI(const FieldsFV3JEDI & other, const oops::Variables & vars)
  : geom_(other.geom_), vars_(vars), time_(other.time_)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  fv3jedi_field_create_f90(keyFlds_, geom_->toFortran(), &conf);
  fv3jedi_field_copy_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsFV3JEDI::~FieldsFV3JEDI() {
  fv3jedi_field_delete_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
FieldsFV3JEDI & FieldsFV3JEDI::operator=(const FieldsFV3JEDI & rhs) {
  fv3jedi_field_copy_f90(keyFlds_, rhs.keyFlds_);
  time_ = rhs.time_;
  return *this;
}
// -----------------------------------------------------------------------------
FieldsFV3JEDI & FieldsFV3JEDI::operator+=(const FieldsFV3JEDI & rhs) {
  fv3jedi_field_self_add_f90(keyFlds_, rhs.keyFlds_);
  return *this;
}
// -----------------------------------------------------------------------------
FieldsFV3JEDI & FieldsFV3JEDI::operator-=(const FieldsFV3JEDI & rhs) {
  fv3jedi_field_self_sub_f90(keyFlds_, rhs.keyFlds_);
  return *this;
}
// -----------------------------------------------------------------------------
FieldsFV3JEDI & FieldsFV3JEDI::operator*=(const double & zz) {
  fv3jedi_field_self_mul_f90(keyFlds_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::zero() {
  fv3jedi_field_zero_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::dirac(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  fv3jedi_field_dirac_f90(keyFlds_, &conf);
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::zero(const util::DateTime & time) {
  fv3jedi_field_zero_f90(keyFlds_);
  time_ = time;
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::axpy(const double & zz, const FieldsFV3JEDI & rhs) {
  fv3jedi_field_axpy_f90(keyFlds_, zz, rhs.keyFlds_);
}
// -----------------------------------------------------------------------------
double FieldsFV3JEDI::dot_product_with(const FieldsFV3JEDI & fld2) const {
  double zz;
  fv3jedi_field_dot_prod_f90(keyFlds_, fld2.keyFlds_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::schur_product_with(const FieldsFV3JEDI & dx) {
    fv3jedi_field_self_schur_f90(keyFlds_, dx.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::random() {
  fv3jedi_field_random_f90(keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::interpolate(const ufo::Locations & locs, const oops::Variables & vars,
                              ufo::GeoVaLs & gom) const {
  const eckit::Configuration * conf = &vars.toFortran();
  fv3jedi_field_interp_f90(keyFlds_, locs.toFortran(), &conf, gom.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::interpolateTL(const ufo::Locations & locs, const oops::Variables & vars,
                                ufo::GeoVaLs & gom) const {
  const eckit::Configuration * conf = &vars.toFortran();
  fv3jedi_field_interp_tl_f90(keyFlds_, locs.toFortran(), &conf, gom.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::interpolateAD(const ufo::Locations & locs, const oops::Variables & vars,
                                const ufo::GeoVaLs & gom) {
  const eckit::Configuration * conf = &vars.toFortran();
  fv3jedi_field_interp_ad_f90(keyFlds_, locs.toFortran(), &conf, gom.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::changeResolution(const FieldsFV3JEDI & other) {
  fv3jedi_field_change_resol_f90(keyFlds_, other.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::add(const FieldsFV3JEDI & rhs) {
  fv3jedi_field_add_incr_f90(keyFlds_, rhs.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::diff(const FieldsFV3JEDI & x1, const FieldsFV3JEDI & x2) {
  fv3jedi_field_diff_incr_f90(keyFlds_, x1.keyFlds_, x2.keyFlds_);
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::convert_to(oops::UnstructuredGrid & ug) const {
  fv3jedi_field_convert_to_f90(keyFlds_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::convert_from(const oops::UnstructuredGrid & ug) {
  fv3jedi_field_convert_from_f90(keyFlds_, ug.toFortran());
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::read(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  fv3jedi_field_read_file_f90(keyFlds_, &conf, &dtp);
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::analytic_init(const eckit::Configuration & config,
				  const GeometryFV3JEDI & geom) {
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  fv3jedi_field_analytic_init_f90(keyFlds_, geom.toFortran(), &conf, &dtp);
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::write(const eckit::Configuration & config) const {
  const eckit::Configuration * conf = &config;
  const util::DateTime * dtp = &time_;
  fv3jedi_field_write_file_f90(keyFlds_, &conf, &dtp);
}
// -----------------------------------------------------------------------------
double FieldsFV3JEDI::norm() const {
  double zz = 0.0;
  fv3jedi_field_rms_f90(keyFlds_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void FieldsFV3JEDI::print(std::ostream & os) const {
  int nx = -1;
  int ny = -1;
  int nf = -1;
  int nb = -1;
  fv3jedi_field_sizes_f90(keyFlds_, nx, ny, nf, nb);
  os << std::endl << "  Resolution = " << nx << ", " << ny
     << ", Fields = " << nf << ", " << nb;
  nf += nb;
  std::vector<double> zstat(3*nf);
  fv3jedi_field_gpnorm_f90(keyFlds_, nf, zstat[0]);
  for (int jj = 0; jj < nf; ++jj) {
    os << std::endl << "  Min=" << zstat[3*jj]
       << ", Max=" << zstat[3*jj+1] << ", RMS=" << zstat[3*jj+2];
  }
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
