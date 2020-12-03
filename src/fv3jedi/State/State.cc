/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------
State::State(const Geometry & geom, const oops::Variables & vars, const util::DateTime & time):
             geom_(new Geometry(geom)), vars_(vars), time_(time) {
  oops::Log::trace() << "State::State (from geom, vars and time) starting" << std::endl;
  fv3jedi_state_create_f90(keyState_, geom_->toFortran(), vars_);
  oops::Log::trace() << "State::State (from geom, vars and time) done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
State::State(const Geometry & geom, const eckit::Configuration & conf): geom_(new Geometry(geom)),
             time_(util::DateTime()) {
  oops::Log::trace() << "State::State (from geom and config) starting" << std::endl;

// Should check if this can be done inside read
  oops::Variables lvars(conf, "state variables");
  this->vars_ = lvars;

  fv3jedi_state_create_f90(keyState_, geom_->toFortran(), vars_);

  // Analytical or read from file
  if (conf.has("analytic_init")) {
    this->analytic_init(conf, geom);
  } else {
    this->read(conf);
  }

  oops::Log::trace() << "State::State (from geom and config) done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
State::State(const Geometry & resol, const State & other): geom_(new Geometry(resol)),
             vars_(other.vars_), time_(other.time_) {
  oops::Log::trace() << "State::State (from geom and other) starting" << std::endl;
  fv3jedi_state_create_f90(keyState_, geom_->toFortran(), vars_);
  fv3jedi_state_change_resol_f90(keyState_, geom_->toFortran(), other.keyState_,
                                 other.geom_->toFortran());
  oops::Log::trace() << "State::State (from geom and other) done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
State::State(const State & other): geom_(other.geom_), vars_(other.vars_), time_(other.time_) {
  oops::Log::trace() << "State::State (from other) starting" << std::endl;
  fv3jedi_state_create_f90(keyState_, geom_->toFortran(), vars_);
  fv3jedi_state_copy_f90(keyState_, other.keyState_);
  oops::Log::trace() << "State::State (from other) done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
State::~State() {
  fv3jedi_state_delete_f90(keyState_);
}
// -------------------------------------------------------------------------------------------------
State & State::operator=(const State & rhs) {
  fv3jedi_state_copy_f90(keyState_, rhs.keyState_);
  time_ = rhs.time_;
  return *this;
}
// -------------------------------------------------------------------------------------------------
void State::changeResolution(const State & other) {
  fv3jedi_state_change_resol_f90(keyState_, geom_->toFortran(), other.keyState_,
                                 other.geom_->toFortran());
}
// -------------------------------------------------------------------------------------------------
State & State::operator+=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());
  // Interpolate increment to state resolution
  Increment dx_sr(*geom_, dx);
  // Call transform and add
  fv3jedi_state_add_incr_f90(geom_->toFortran(), keyState_, dx_sr.toFortran());
  return *this;
}
// -------------------------------------------------------------------------------------------------
void State::read(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  fv3jedi_state_read_file_f90(geom_->toFortran(), keyState_, &conf, &dtp);
}
// -------------------------------------------------------------------------------------------------
void State::analytic_init(const eckit::Configuration & config, const Geometry & geom) {
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  fv3jedi_state_analytic_init_f90(keyState_, geom.toFortran(), &conf, &dtp);
}
// -------------------------------------------------------------------------------------------------
void State::write(const eckit::Configuration & config) const {
  const eckit::Configuration * conf = &config;
  const util::DateTime * dtp = &time_;
  fv3jedi_state_write_file_f90(geom_->toFortran(), keyState_, &conf, &dtp);
}
// -------------------------------------------------------------------------------------------------
void State::print(std::ostream & os) const {
  // Get the number of fields
  int numberFields;
  int cubeSize;
  fv3jedi_state_getnfieldsncube_f90(keyState_, numberFields, cubeSize);

  // Header
  os << std::endl
     << " -----------------------------------------------"
        "------------------------------------------------";
  os << std::endl << " State print | number of fields = " << numberFields
                  << " | cube sphere face size: C" << cubeSize;

  // Print info field by field
  const int FieldNameLen = 31;
  char fieldName[FieldNameLen];
  std::vector<double> minMaxRms(3);
  for (int f = 0; f < numberFields; f++) {
    int fp1 = f+1;
    fv3jedi_state_getminmaxrms_f90(keyState_, fp1, FieldNameLen, fieldName, minMaxRms[0]);
    std::string fieldNameStr(fieldName);
    os << std::endl << std::scientific << std::showpos << "   "
                    << fieldNameStr.substr(0, FieldNameLen-1) << ": Min = " << minMaxRms[0]
                    << ", Max = " << minMaxRms[1] << ", RMS = " << minMaxRms[2]
                    << std::noshowpos;  //  << std::defaultfloat;
  }

  os.unsetf(std::ios_base::floatfield);

  // Footer
  os << std::endl
     << " -----------------------------------------------"
        "------------------------------------------------";
}
// -----------------------------------------------------------------------------
void State::zero() {
  fv3jedi_state_zero_f90(keyState_);
}
// -----------------------------------------------------------------------------
void State::accumul(const double & zz, const State & xx) {
  fv3jedi_state_axpy_f90(keyState_, zz, xx.keyState_);
}
// -----------------------------------------------------------------------------
double State::norm() const {
  double zz = 0.0;
  fv3jedi_state_norm_f90(keyState_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
size_t State::serialSize() const {
  oops::Log::trace() << "State serialSize starting" << std::endl;
  size_t nn = 1;
  int sz = 0;
  fv3jedi_state_sersize_f90(keyState_, sz);
  nn += sz;
  nn += time_.serialSize();
  return nn;
  oops::Log::trace() << "State serialSize done" << std::endl;
}
// -----------------------------------------------------------------------------
void State::serialize(std::vector<double> & vect) const {
  oops::Log::trace() << "State serialize starting" << std::endl;
  int size_fld = this->serialSize() - 3;
  std::vector<double> v_fld(size_fld, 0);

  fv3jedi_state_serialize_f90(keyState_, size_fld, v_fld.data());
  vect.insert(vect.end(), v_fld.begin(), v_fld.end());

  // Serialize the date and time
  vect.push_back(-54321.56789);
  time_.serialize(vect);

  oops::Log::trace() << "State serialize done" << std::endl;
}
// -----------------------------------------------------------------------------
void State::deserialize(const std::vector<double> & vect, size_t & index) {
  oops::Log::trace() << "State deserialize starting" << std::endl;
  fv3jedi_state_deserialize_f90(keyState_, vect.size(), vect.data(), index);

  ASSERT(vect.at(index) == -54321.56789);
  ++index;

  time_.deserialize(vect, index);
  oops::Log::trace() << "State deserialize done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
