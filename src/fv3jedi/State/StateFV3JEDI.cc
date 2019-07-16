/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <algorithm>
#include <string>
#include <vector>


#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

#include "fv3jedi/Geometry/GeometryFV3JEDI.h"
#include "fv3jedi/GetValues/GetValuesTrajFV3JEDI.h"
#include "fv3jedi/Increment/IncrementFV3JEDI.h"
#include "fv3jedi/State/StateFV3JEDI.h"
#include "fv3jedi/Utilities/UtilitiesFV3JEDI.h"
#include "StateFV3JEDIFortran.h"

namespace fv3jedi {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
StateFV3JEDI::StateFV3JEDI(const GeometryFV3JEDI & geom,
                           const oops::Variables & vars,
                           const util::DateTime & time):
  geom_(new GeometryFV3JEDI(geom)), vars_(vars), time_(time)
{
  const eckit::Configuration * confvars = &vars_.toFortran();
  fv3jedi_state_create_f90(keyState_, geom_->toFortran(), &confvars);
  oops::Log::trace() << "StateFV3JEDI::StateFV3JEDI created." << std::endl;
}
// -----------------------------------------------------------------------------
StateFV3JEDI::StateFV3JEDI(const GeometryFV3JEDI & geom,
                           const oops::Variables & vars,
                           const eckit::Configuration & sconf)
  : geom_(new GeometryFV3JEDI(geom)), time_(util::DateTime())
{
  oops::Log::trace() << "StateFV3JEDI::StateFV3JEDI create from analytical or"
                        " from file." << std::endl;
  // Optionally user can overwrite incoming variables from config
  oops::Variables lvars;
  if (sconf.has("variables")) {
    oops::Variables lvars(sconf);
    this->vars_ = lvars;
  } else {
    this->vars_ = vars;
  }
  const eckit::Configuration * confvars = &vars_.toFortran();
  fv3jedi_state_create_f90(keyState_, geom_->toFortran(), &confvars);

  const eckit::Configuration * conf = &sconf;
  util::DateTime * dtp = &time_;
  if (sconf.has("analytic_init")) {
    stageFv3Files(sconf);
    fv3jedi_state_analytic_init_f90(keyState_, geom.toFortran(), &conf, &dtp);
    removeFv3Files();
  } else {
    fv3jedi_state_read_file_f90(geom_->toFortran(), keyState_, &conf, &dtp);
  }

  oops::Log::trace() << "StateFV3JEDI::StateFV3JEDI create from analytical or"
                        " from file done." << std::endl;
}
// -----------------------------------------------------------------------------
StateFV3JEDI::StateFV3JEDI(const GeometryFV3JEDI & resol,
                           const StateFV3JEDI & other):
  geom_(new GeometryFV3JEDI(resol)), vars_(other.vars_), time_(other.time_)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  fv3jedi_state_create_f90(keyState_, geom_->toFortran(), &conf);
  fv3jedi_state_change_resol_f90(keyState_, geom_->toFortran(), other.keyState_,
                                 other.geom_->toFortran());
  oops::Log::trace() << "StateFV3JEDI::StateFV3JEDI created by interpolation."
                     << std::endl;
}
// -----------------------------------------------------------------------------
StateFV3JEDI::StateFV3JEDI(const StateFV3JEDI & other):
  geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  const eckit::Configuration * conf = &vars_.toFortran();
  fv3jedi_state_create_f90(keyState_, geom_->toFortran(), &conf);
  fv3jedi_state_copy_f90(keyState_, other.keyState_);
  oops::Log::trace() << "StateFV3JEDI::StateFV3JEDI copied." << std::endl;
}
// -----------------------------------------------------------------------------
StateFV3JEDI::~StateFV3JEDI() {
  fv3jedi_state_delete_f90(keyState_);
  oops::Log::trace() << "StateFV3JEDI::StateFV3JEDI destructed." << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
StateFV3JEDI & StateFV3JEDI::operator=(const StateFV3JEDI & rhs) {
  fv3jedi_state_copy_f90(keyState_, rhs.keyState_);
  time_ = rhs.time_;
  return *this;
}
// -----------------------------------------------------------------------------
/// Get state values at observation locations
// -----------------------------------------------------------------------------
void StateFV3JEDI::getValues(const ufo::Locations & locs,
                             const oops::Variables & vars,
                             ufo::GeoVaLs & gom) const {
  oops::Log::trace() << "StateFV3JEDI::getValues starting." << std::endl;
  const eckit::Configuration * conf = &vars.toFortran();
  fv3jedi_state_getvalues_notraj_f90(geom_->toFortran(), keyState_,
                                     locs.toFortran(), &conf,
                                     gom.toFortran());
  oops::Log::trace() << "StateFV3JEDI::getValues done." << std::endl;
}
// -----------------------------------------------------------------------------
void StateFV3JEDI::getValues(const ufo::Locations & locs,
                             const oops::Variables & vars,
                             ufo::GeoVaLs & gom,
                             const GetValuesTrajFV3JEDI & traj) const {
  oops::Log::trace() << "StateFV3JEDI::getValues traj starting." << std::endl;
  const eckit::Configuration * conf = &vars.toFortran();
  fv3jedi_state_getvalues_f90(geom_->toFortran(), keyState_, locs.toFortran(),
                              &conf, gom.toFortran(), traj.toFortran());
  oops::Log::trace() << "StateFV3JEDI::getValues traj done." << std::endl;
}
// -----------------------------------------------------------------------------
/// Interpolate full state
// -----------------------------------------------------------------------------
void StateFV3JEDI::changeResolution(const StateFV3JEDI & other) {
  oops::Log::trace() << "StateFV3JEDI change resolution starting" << std::endl;
  fv3jedi_state_change_resol_f90(keyState_, geom_->toFortran(), other.keyState_,
                                 other.geom_->toFortran());
  oops::Log::trace() << "StateFV3JEDI change resolution done" << std::endl;
}
// -----------------------------------------------------------------------------
/// Interactions with Increments
// -----------------------------------------------------------------------------
StateFV3JEDI & StateFV3JEDI::operator+=(const IncrementFV3JEDI & dx) {
  oops::Log::trace() << "StateFV3JEDI add increment starting" << std::endl;
  ASSERT(this->validTime() == dx.validTime());
  // Interpolate increment to state resolution
  IncrementFV3JEDI dx_sr(*geom_, dx);
  // Call transform and add
  fv3jedi_state_add_incr_f90(geom_->toFortran(), keyState_, dx_sr.toFortran());
  oops::Log::trace() << "StateFV3JEDI add increment done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void StateFV3JEDI::read(const eckit::Configuration & config) {
  oops::Log::trace() << "StateFV3JEDI read starting" << std::endl;
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  fv3jedi_state_read_file_f90(geom_->toFortran(), keyState_, &conf, &dtp);
  oops::Log::trace() << "StateFV3JEDI read done" << std::endl;
}
// -----------------------------------------------------------------------------
void StateFV3JEDI::analytic_init(const eckit::Configuration & config,
                                 const GeometryFV3JEDI & geom) {
  oops::Log::trace() << "StateFV3JEDI analytic init starting" << std::endl;
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  stageFv3Files(config);
  fv3jedi_state_analytic_init_f90(keyState_, geom.toFortran(), &conf, &dtp);
  removeFv3Files();
  oops::Log::trace() << "StateFV3JEDI analytic init done" << std::endl;
}
// -----------------------------------------------------------------------------
void StateFV3JEDI::write(const eckit::Configuration & config) const {
  oops::Log::trace() << "StateFV3JEDI write starting" << std::endl;
  const eckit::Configuration * conf = &config;
  const util::DateTime * dtp = &time_;
  fv3jedi_state_write_file_f90(geom_->toFortran(), keyState_, &conf, &dtp);
  oops::Log::trace() << "StateFV3JEDI write done" << std::endl;
}
// -----------------------------------------------------------------------------
void StateFV3JEDI::print(std::ostream & os) const {
  oops::Log::trace() << "StateFV3JEDI print starting" << std::endl;
  fv3jedi_state_print_f90(keyState_);
//  oops::Log::debug() << "  Valid time: " << validTime() << std::endl;
//  int nx = 0;
//  int ny = 0;
//  int nf = 8;
//  fv3jedi_state_sizes_f90(keyState_, nx, ny, nf);
//  oops::Log::debug() << "Cube faces = " << nx << "x" << ny
//     << ", Number of state fields = " << nf << std::endl;
//  std::vector<double> zstat(3*nf);
//  fv3jedi_state_gpnorm_f90(keyState_, nf, zstat[0]);
//  for (int jj = 0; jj < nf; ++jj) {
//    oops::Log::debug() << "State=" << jj+1 <<"  Min=" << zstat[3*jj]
//       << ", Max=" << zstat[3*jj+1] << ", RMS=" << zstat[3*jj+2] << std::endl;
//  }
  oops::Log::trace() << "StateFV3JEDI print done" << std::endl;
}
// -----------------------------------------------------------------------------
/// For accumulator
// -----------------------------------------------------------------------------
void StateFV3JEDI::zero() {
  oops::Log::trace() << "StateFV3JEDI zero starting" << std::endl;
  fv3jedi_state_zero_f90(keyState_);
  oops::Log::trace() << "StateFV3JEDI zero done" << std::endl;
}
// -----------------------------------------------------------------------------
void StateFV3JEDI::accumul(const double & zz, const StateFV3JEDI & xx) {
  oops::Log::trace() << "StateFV3JEDI accumul starting" << std::endl;
  fv3jedi_state_axpy_f90(keyState_, zz, xx.keyState_);
  oops::Log::trace() << "StateFV3JEDI accumul done" << std::endl;
}
// -----------------------------------------------------------------------------
double StateFV3JEDI::norm() const {
  oops::Log::trace() << "StateFV3JEDI norm starting" << std::endl;
  double zz = 0.0;
  fv3jedi_state_rms_f90(keyState_, zz);
  return zz;
  oops::Log::trace() << "StateFV3JEDI norm done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
