/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <iomanip>
#include <ios>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/Increment/Increment.h"
#include "fv3jedi/State/State.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------
Increment::Increment(const Geometry & geom, const oops::Variables & vars,
                     const util::DateTime & time): geom_(new Geometry(geom)), vars_(vars),
                     time_(time) {
  oops::Log::trace() << "Increment::Increment (from geom, vars and time) starting" << std::endl;
  fv3jedi_increment_create_f90(keyInc_, geom_->toFortran(), vars_);
  fv3jedi_increment_zero_f90(keyInc_);
  oops::Log::trace() << "Increment::Increment (from geom, vars and time) done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
Increment::Increment(const Geometry & geom, const Increment & other): geom_(new Geometry(geom)),
                     vars_(other.vars_), time_(other.time_) {
  oops::Log::trace() << "Increment::Increment (from geom and other) starting" << std::endl;
  fv3jedi_increment_create_f90(keyInc_, geom_->toFortran(), vars_);
  fv3jedi_increment_change_resol_f90(keyInc_, geom_->toFortran(), other.keyInc_,
                                     other.geometry()->toFortran());
  oops::Log::trace() << "Increment::Increment (from geom and other) done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
Increment::Increment(const Increment & other, const bool copy): geom_(other.geom_),
                     vars_(other.vars_), time_(other.time_) {
  oops::Log::trace() << "Increment::Increment (from other and bool copy) starting" << std::endl;
  fv3jedi_increment_create_f90(keyInc_, geom_->toFortran(), vars_);
  if (copy) {
    fv3jedi_increment_copy_f90(keyInc_, other.keyInc_);
  } else {
    fv3jedi_increment_zero_f90(keyInc_);
  }
  oops::Log::trace() << "Increment::Increment (from other and bool copy) done" << std::endl;
}
// -------------------------------------------------------------------------------------------------
Increment::~Increment() {
  fv3jedi_increment_delete_f90(keyInc_);
}
// -------------------------------------------------------------------------------------------------
void Increment::diff(const State & x1, const State & x2) {
  ASSERT(this->validTime() == x1.validTime());
  ASSERT(this->validTime() == x2.validTime());
  // States at increment resolution
  State x1_ir(*geom_, x1);
  State x2_ir(*geom_, x2);
  fv3jedi_increment_diff_incr_f90(keyInc_, x1_ir.toFortran(), x2_ir.toFortran(),
                                  geom_->toFortran());
}
// -------------------------------------------------------------------------------------------------
Increment & Increment::operator=(const Increment & rhs) {
  fv3jedi_increment_copy_f90(keyInc_, rhs.keyInc_);
  time_ = rhs.time_;
  return *this;
}
// -------------------------------------------------------------------------------------------------
Increment & Increment::operator+=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());
  fv3jedi_increment_self_add_f90(keyInc_, dx.keyInc_);
  return *this;
}
// -------------------------------------------------------------------------------------------------
Increment & Increment::operator-=(const Increment & dx) {
  ASSERT(this->validTime() == dx.validTime());
  fv3jedi_increment_self_sub_f90(keyInc_, dx.keyInc_);
  return *this;
}
// -------------------------------------------------------------------------------------------------
Increment & Increment::operator*=(const double & zz) {
  fv3jedi_increment_self_mul_f90(keyInc_, zz);
  return *this;
}
// -------------------------------------------------------------------------------------------------
void Increment::zero() {
  fv3jedi_increment_zero_f90(keyInc_);
}
// -------------------------------------------------------------------------------------------------
void Increment::zero(const util::DateTime & vt) {
  fv3jedi_increment_zero_f90(keyInc_);
  time_ = vt;
}
// -------------------------------------------------------------------------------------------------
void Increment::ones() {
  fv3jedi_increment_ones_f90(keyInc_);
}
// -------------------------------------------------------------------------------------------------
void Increment::axpy(const double & zz, const Increment & dx, const bool check) {
  ASSERT(!check || this->validTime() == dx.validTime());
  fv3jedi_increment_axpy_inc_f90(keyInc_, zz, dx.keyInc_);
}
// -------------------------------------------------------------------------------------------------
void Increment::accumul(const double & zz, const State & xx) {
  fv3jedi_increment_axpy_state_f90(keyInc_, zz, xx.toFortran());
}
// -------------------------------------------------------------------------------------------------
void Increment::schur_product_with(const Increment & dx) {
  fv3jedi_increment_self_schur_f90(keyInc_, dx.keyInc_);
}
// -------------------------------------------------------------------------------------------------
double Increment::dot_product_with(const Increment & other) const {
  double zz;
  fv3jedi_increment_dot_prod_f90(keyInc_, other.keyInc_, zz);
  return zz;
}
// -------------------------------------------------------------------------------------------------
void Increment::random() {
  fv3jedi_increment_random_f90(keyInc_);
}
// -------------------------------------------------------------------------------------------------
oops::LocalIncrement Increment::getLocal(const GeometryIterator & iter) const {
  int ist, iend, jst, jend, npz;
  fv3jedi_geom_start_end_f90(geom_->toFortran(), ist, iend, jst, jend, npz);

  oops::Variables fieldNames = vars_;

  std::vector<int> varlens(fieldNames.size());
  for (unsigned int ii = 0; ii < fieldNames.size(); ii++) {
     // might need to modify if non-npz variables are also used
     varlens[ii] = npz;
  }

  int lenvalues = std::accumulate(varlens.begin(), varlens.end(), 0);
  std::vector<double> values(lenvalues);

  // Get variable values
  fv3jedi_increment_getpoint_f90(keyInc_, iter.toFortran(), values[0],
                                  values.size());

  return oops::LocalIncrement(oops::Variables(fieldNames), values, varlens);
}
// -------------------------------------------------------------------------------------------------
void Increment::setLocal(const oops::LocalIncrement & values, const GeometryIterator & iter) {
  const std::vector<double> vals = values.getVals();
  fv3jedi_increment_setpoint_f90(keyInc_, iter.toFortran(), vals[0], vals.size());
}
// -------------------------------------------------------------------------------------------------
void Increment::setAtlas(atlas::FieldSet * afieldset) const {
  fv3jedi_increment_set_atlas_f90(keyInc_, geom_->toFortran(), vars_, afieldset->get());
}
// -------------------------------------------------------------------------------------------------
void Increment::toAtlas(atlas::FieldSet * afieldset) const {
  fv3jedi_increment_to_atlas_f90(keyInc_, geom_->toFortran(), vars_, afieldset->get());
}
// -------------------------------------------------------------------------------------------------
void Increment::fromAtlas(atlas::FieldSet * afieldset) {
  fv3jedi_increment_from_atlas_f90(keyInc_, geom_->toFortran(), vars_, afieldset->get());
}
// -------------------------------------------------------------------------------------------------
void Increment::read(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  util::DateTime * dtp = &time_;
  fv3jedi_increment_read_file_f90(geom_->toFortran(), keyInc_, &conf, &dtp);
}
// -------------------------------------------------------------------------------------------------
void Increment::write(const eckit::Configuration & config) const {
  const eckit::Configuration * conf = &config;
  const util::DateTime * dtp = &time_;
  fv3jedi_increment_write_file_f90(geom_->toFortran(), keyInc_, &conf, &dtp);
}
// -------------------------------------------------------------------------------------------------
double Increment::norm() const {
  double zz = 0.0;
  fv3jedi_increment_norm_f90(keyInc_, zz);
  return zz;
}
// -------------------------------------------------------------------------------------------------
void Increment::print(std::ostream & os) const {
  // Get the number of fields
  int numberFields;
  int cubeSize;
  fv3jedi_increment_getnfieldsncube_f90(keyInc_, numberFields, cubeSize);

  // Header
  os << std::endl
     << " -----------------------------------------------"
        "------------------------------------------------";
  os << std::endl << " Increment print | number of fields = " << numberFields
                  << " | cube sphere face size: C" << cubeSize;

  // Print info field by field
  const int FieldNameLen = 31;
  char fieldName[FieldNameLen];
  std::vector<double> minMaxRms(3);
  for (int f = 0; f < numberFields; f++) {
    int fp1 = f+1;
    fv3jedi_increment_getminmaxrms_f90(keyInc_, fp1, FieldNameLen, fieldName, minMaxRms[0]);
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
// -------------------------------------------------------------------------------------------------
void Increment::dirac(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  fv3jedi_increment_dirac_f90(keyInc_, &conf, geom_->toFortran());
}
// -------------------------------------------------------------------------------------------------
size_t Increment::serialSize() const {
  size_t nn = 1;
  int inc_size;
  fv3jedi_increment_sizes_f90(keyInc_, inc_size);
  nn+= inc_size;  // to verify
  nn += time_.serialSize();
  return nn;
}
// -------------------------------------------------------------------------------------------------
void Increment::serialize(std::vector<double> & vect) const {
  int size_inc = this->serialSize() - 3;
  std::vector<double> v_inc(size_inc, 0);

  fv3jedi_increment_serialize_f90(keyInc_, size_inc, v_inc.data());
  vect.insert(vect.end(), v_inc.begin(), v_inc.end());

  // Serialize the date and time
  vect.push_back(-54321.98765);
  time_.serialize(vect);
}
// -------------------------------------------------------------------------------------------------
void Increment::deserialize(const std::vector<double> & vect,
                                   size_t & index) {
  fv3jedi_increment_deserialize_f90(keyInc_, vect.size(), vect.data(), index);

  ASSERT(vect.at(index) == -54321.98765);
  ++index;

  time_.deserialize(vect, index);
}
// -------------------------------------------------------------------------------------------------
}  // namespace fv3jedi
