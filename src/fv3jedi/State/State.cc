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
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/GetValues/GetValuesTraj.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"
#include "fv3jedi/Utilities/Utilities.h"

namespace fv3jedi {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
State::State(const Geometry & geom,
                           const oops::Variables & vars,
                           const util::DateTime & time):
  geom_(new Geometry(geom)), vars_(vars), time_(time)
{
  fv3jedi_state_create_f90(keyState_, geom_->toFortran(), vars_);
  oops::Log::trace() << "State::State created." << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const Geometry & geom,
                           const oops::Variables & vars,
                           const eckit::Configuration & sconf)
  : geom_(new Geometry(geom)), time_(util::DateTime())
{
  oops::Log::trace() << "State::State create from analytical or"
                        " from file." << std::endl;
  // Optionally user can overwrite incoming variables from config
  oops::Variables lvars;
  if (sconf.has("variables")) {
    oops::Variables lvars(sconf);
    this->vars_ = lvars;
  } else {
    this->vars_ = vars;
  }
  fv3jedi_state_create_f90(keyState_, geom_->toFortran(), vars_);

  const eckit::Configuration * conf = &sconf;
  util::DateTime * dtp = &time_;
  if (sconf.has("analytic_init")) {
    stageFv3Files(sconf);
    fv3jedi_state_analytic_init_f90(keyState_, geom.toFortran(), &conf, &dtp);
    removeFv3Files();
  } else {
    fv3jedi_state_read_file_f90(geom_->toFortran(), keyState_, &conf, &dtp);
  }

  oops::Log::trace() << "State::State create from analytical or"
                        " from file done." << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const Geometry & resol,
                           const State & other):
  geom_(new Geometry(resol)), vars_(other.vars_), time_(other.time_)
{
  fv3jedi_state_create_f90(keyState_, geom_->toFortran(), vars_);
  fv3jedi_state_change_resol_f90(keyState_, geom_->toFortran(), other.keyState_,
                                 other.geom_->toFortran());
  oops::Log::trace() << "State::State created by interpolation."
                     << std::endl;
}
// -----------------------------------------------------------------------------
State::State(const State & other):
  geom_(other.geom_), vars_(other.vars_), time_(other.time_)
{
  fv3jedi_state_create_f90(keyState_, geom_->toFortran(), vars_);
  fv3jedi_state_copy_f90(keyState_, other.keyState_);
  oops::Log::trace() << "State::State copied." << std::endl;
}
// -----------------------------------------------------------------------------
State::~State() {
  fv3jedi_state_delete_f90(keyState_);
  oops::Log::trace() << "State::State destructed." << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
State & State::operator=(const State & rhs) {
  fv3jedi_state_copy_f90(keyState_, rhs.keyState_);
  time_ = rhs.time_;
  return *this;
}
// -----------------------------------------------------------------------------
/// Get state values at observation locations
// -----------------------------------------------------------------------------
void State::getValues(const ufo::Locations & locs,
                             const oops::Variables & vars,
                             ufo::GeoVaLs & gom) const {
  oops::Log::trace() << "State::getValues starting." << std::endl;
  fv3jedi_state_getvalues_notraj_f90(geom_->toFortran(), keyState_,
                                     locs.toFortran(), vars,
                                     gom.toFortran());
  oops::Log::trace() << "State::getValues done." << std::endl;
}
// -----------------------------------------------------------------------------
void State::getValues(const ufo::Locations & locs,
                             const oops::Variables & vars,
                             ufo::GeoVaLs & gom,
                             const GetValuesTrajMatrix & traj) const {
  oops::Log::trace() << "State::getValues traj starting." << std::endl;
  fv3jedi_state_getvalues_f90(geom_->toFortran(), keyState_, locs.toFortran(),
                              vars, gom.toFortran(), traj.toFortran());
  oops::Log::trace() << "State::getValues traj done." << std::endl;
}
// -----------------------------------------------------------------------------
/// Interpolate full state
// -----------------------------------------------------------------------------
void State::changeResolution(const State & other) {
  oops::Log::trace() << "State change resolution starting" << std::endl;
  fv3jedi_state_change_resol_f90(keyState_, geom_->toFortran(), other.keyState_,
                                 other.geom_->toFortran());
  oops::Log::trace() << "State change resolution done" << std::endl;
}
// -----------------------------------------------------------------------------
/// Interactions with Increments
// -----------------------------------------------------------------------------
State & State::operator+=(const Increment & dx) {
  oops::Log::trace() << "State add increment starting" << std::endl;
  ASSERT(this->validTime() == dx.validTime());
  // Interpolate increment to state resolution
  Increment dx_sr(*geom_, dx);
  // Call transform and add
  fv3jedi_state_add_incr_f90(geom_->toFortran(), keyState_, dx_sr.toFortran());
  oops::Log::trace() << "State add increment done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void State::read(const eckit::Configuration & config) {
  oops::Log::trace() << "State read starting" << std::endl;
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  fv3jedi_state_read_file_f90(geom_->toFortran(), keyState_, &conf, &dtp);
  oops::Log::trace() << "State read done" << std::endl;
}
// -----------------------------------------------------------------------------
void State::analytic_init(const eckit::Configuration & config,
                                 const Geometry & geom) {
  oops::Log::trace() << "State analytic init starting" << std::endl;
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  stageFv3Files(config);
  fv3jedi_state_analytic_init_f90(keyState_, geom.toFortran(), &conf, &dtp);
  removeFv3Files();
  oops::Log::trace() << "State analytic init done" << std::endl;
}
// -----------------------------------------------------------------------------
void State::write(const eckit::Configuration & config) const {
  oops::Log::trace() << "State write starting" << std::endl;
  const eckit::Configuration * conf = &config;
  const util::DateTime * dtp = &time_;
  fv3jedi_state_write_file_f90(geom_->toFortran(), keyState_, &conf, &dtp);
  oops::Log::trace() << "State write done" << std::endl;
}
// -----------------------------------------------------------------------------
void State::print(std::ostream & os) const {
  oops::Log::trace() << "State print starting" << std::endl;
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
  oops::Log::trace() << "State print done" << std::endl;
}
// -----------------------------------------------------------------------------
/// For accumulator
// -----------------------------------------------------------------------------
void State::zero() {
  oops::Log::trace() << "State zero starting" << std::endl;
  fv3jedi_state_zero_f90(keyState_);
  oops::Log::trace() << "State zero done" << std::endl;
}
// -----------------------------------------------------------------------------
void State::accumul(const double & zz, const State & xx) {
  oops::Log::trace() << "State accumul starting" << std::endl;
  fv3jedi_state_axpy_f90(keyState_, zz, xx.keyState_);
  oops::Log::trace() << "State accumul done" << std::endl;
}
// -----------------------------------------------------------------------------
double State::norm() const {
  oops::Log::trace() << "State norm starting" << std::endl;
  double zz = 0.0;
  fv3jedi_state_rms_f90(keyState_, zz);
  return zz;
  oops::Log::trace() << "State norm done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
