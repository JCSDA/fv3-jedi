/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <algorithm>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ioda/Locations.h"
#include "ufo/GeoVaLs.h"

#include "src/Increment/IncrementFV3JEDI.h"
#include "ErrorCovarianceFV3JEDI.h"
#include "GeometryFV3JEDI.h"
#include "GetValuesTrajFV3JEDI.h"
#include "IncrementFV3JEDIFortran.h"
#include "StateFV3JEDI.h"

namespace fv3jedi {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
IncrementFV3JEDI::IncrementFV3JEDI(const GeometryFV3JEDI & geom,
                                   const oops::Variables & vars,
                                   const util::DateTime & time):
  geom_(new GeometryFV3JEDI(geom)), vars_(vars), time_(time)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  fv3jedi_increment_create_f90(keyInc_, geom_->toFortran(), &conf);
  fv3jedi_increment_zero_f90(keyInc_);
  oops::Log::trace() << "IncrementFV3JEDI constructed." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementFV3JEDI::IncrementFV3JEDI(const GeometryFV3JEDI & geom,
                                   const IncrementFV3JEDI & other)
  : geom_(new GeometryFV3JEDI(geom)), vars_(other.vars_), time_(other.time_)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  fv3jedi_increment_create_f90(keyInc_, geom_->toFortran(), &conf);
  fv3jedi_increment_change_resol_f90(keyInc_, other.keyInc_);
  oops::Log::trace() << "IncrementFV3JEDI constructed from other." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementFV3JEDI::IncrementFV3JEDI(const IncrementFV3JEDI & other,
                                   const bool copy)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  fv3jedi_increment_create_f90(keyInc_, geom_->toFortran(), &conf);
  if (copy) {
    fv3jedi_increment_copy_f90(keyInc_, other.keyInc_);
  } else {
    fv3jedi_increment_zero_f90(keyInc_);
  }
  oops::Log::trace() << "IncrementFV3JEDI copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementFV3JEDI::IncrementFV3JEDI(const IncrementFV3JEDI & other)
  : geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  fv3jedi_increment_create_f90(keyInc_, geom_->toFortran(), &conf);
  fv3jedi_increment_copy_f90(keyInc_, other.keyInc_);
  oops::Log::trace() << "IncrementFV3JEDI copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementFV3JEDI::~IncrementFV3JEDI() {
  fv3jedi_increment_delete_f90(keyInc_);
  oops::Log::trace() << "IncrementFV3JEDI destructed" << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::diff(const StateFV3JEDI & x1, const StateFV3JEDI & x2) {
  ASSERT(this->validTime() == x1.validTime());
  ASSERT(this->validTime() == x2.validTime());
  fv3jedi_increment_diff_incr_f90(keyInc_, x1.toFortran(), x2.toFortran());
}
// -----------------------------------------------------------------------------
IncrementFV3JEDI & IncrementFV3JEDI::operator=(const IncrementFV3JEDI & rhs) {
  fv3jedi_increment_copy_f90(keyInc_, rhs.keyInc_);
  time_ = rhs.time_;
  return *this;
}
// -----------------------------------------------------------------------------
IncrementFV3JEDI & IncrementFV3JEDI::operator+=(const IncrementFV3JEDI & dx) {
  ASSERT(this->validTime() == dx.validTime());
  fv3jedi_increment_self_add_f90(keyInc_, dx.keyInc_);
  return *this;
}
// -----------------------------------------------------------------------------
IncrementFV3JEDI & IncrementFV3JEDI::operator-=(const IncrementFV3JEDI & dx) {
  ASSERT(this->validTime() == dx.validTime());
  fv3jedi_increment_self_sub_f90(keyInc_, dx.keyInc_);
  return *this;
}
// -----------------------------------------------------------------------------
IncrementFV3JEDI & IncrementFV3JEDI::operator*=(const double & zz) {
  fv3jedi_increment_self_mul_f90(keyInc_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::zero() {
  fv3jedi_increment_zero_f90(keyInc_);
}
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::zero(const util::DateTime & vt) {
  fv3jedi_increment_zero_f90(keyInc_);
  time_ = vt;
}
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::axpy(const double & zz, const IncrementFV3JEDI & dx,
                       const bool check) {
  ASSERT(!check || this->validTime() == dx.validTime());
  fv3jedi_increment_axpy_inc_f90(keyInc_, zz, dx.keyInc_);
}
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::accumul(const double & zz, const StateFV3JEDI & xx) {
  fv3jedi_increment_axpy_state_f90(keyInc_, zz, xx.toFortran());
}
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::schur_product_with(const IncrementFV3JEDI & dx) {
  fv3jedi_increment_self_schur_f90(keyInc_, dx.keyInc_);
}
// -----------------------------------------------------------------------------
double IncrementFV3JEDI::dot_product_with(const IncrementFV3JEDI & other)
                                          const {
  double zz;
  fv3jedi_increment_dot_prod_f90(keyInc_, other.keyInc_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::random() {
  fv3jedi_increment_random_f90(keyInc_);
}
// -----------------------------------------------------------------------------
/// Get increment values at observation locations
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::getValuesTL(const ioda::Locations & locs,
                                   const oops::Variables & vars,
                                   ufo::GeoVaLs & gom,
                                   const GetValuesTrajFV3JEDI & traj) const {
  const eckit::Configuration * conf = &vars.toFortran();
  fv3jedi_increment_getvalues_tl_f90(geom_->toFortran(),
                                     keyInc_,
                                     locs.toFortran(),
                                     &conf,
                                     gom.toFortran(),
                                     traj.toFortran());
}
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::getValuesAD(const ioda::Locations & locs,
                             const oops::Variables & vars,
                             const ufo::GeoVaLs & gom,
                             const GetValuesTrajFV3JEDI & traj) {
  const eckit::Configuration * conf = &vars.toFortran();
  fv3jedi_increment_getvalues_ad_f90(geom_->toFortran(),
                                     keyInc_, locs.toFortran(),
                                     &conf,
                                     gom.toFortran(),
                                     traj.toFortran());
}
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::ug_coord(oops::UnstructuredGrid & ug,
                                const int & colocated) const {
  fv3jedi_increment_ug_coord_f90(keyInc_, ug.toFortran(), colocated,
                             geom_->toFortran());
}
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::field_to_ug(oops::UnstructuredGrid & ug,
                                   const int & colocated) const {
  fv3jedi_increment_increment_to_ug_f90(keyInc_, ug.toFortran(), colocated);
}
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::field_from_ug(const oops::UnstructuredGrid & ug) {
  fv3jedi_increment_increment_from_ug_f90(keyInc_, ug.toFortran());
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::read(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  fv3jedi_increment_read_file_f90(geom_->toFortran(), keyInc_, &conf, &dtp);
}
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::write(const eckit::Configuration & config) const {
  const eckit::Configuration * conf = &config;
  const util::DateTime * dtp = &time_;
  fv3jedi_increment_write_file_f90(geom_->toFortran(), keyInc_, &conf, &dtp);
}
// -----------------------------------------------------------------------------
double IncrementFV3JEDI::norm() const {
  double zz = 0.0;
  fv3jedi_increment_rms_f90(keyInc_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::print(std::ostream & os) const {
  oops::Log::trace() << "IncrementFV3JEDI print starting" << std::endl;
  os << std::endl << "  Valid time: " << validTime();
  int nx = 0;
  int ny = 0;
  int nf = 5;
  fv3jedi_increment_sizes_f90(keyInc_, nx, ny, nf);
  os << std::endl << "  Resolution = " << nx << ", " << ny
     << ", Increment = " << nf;
  std::vector<double> zstat(3*nf);
  fv3jedi_increment_gpnorm_f90(keyInc_, nf, zstat[0]);
  for (int jj = 0; jj < nf; ++jj) {
    os << std::endl <<"Increment=" << jj+1 <<"  Min=" << zstat[3*jj]
       << ", Max=" << zstat[3*jj+1] << ", RMS=" << zstat[3*jj+2];
  }
  oops::Log::trace() << "IncrementFV3JEDI print done" << std::endl;
}
// -----------------------------------------------------------------------------
void IncrementFV3JEDI::dirac(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  fv3jedi_increment_dirac_f90(keyInc_, &conf, geom_->toFortran());
}
// -----------------------------------------------------------------------------
}  // namespace fv3jedi
