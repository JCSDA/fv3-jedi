/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <algorithm>
#include <string>
#include <vector>
#include "src/StateFV3JEDI.h"
#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ioda/Locations.h"
#include "ModelBiasFV3JEDI.h"
#include "FieldsFV3JEDI.h"
#include "GeometryFV3JEDI.h"
#include "IncrementFV3JEDI.h"
#include "ModelFV3JEDI.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "GetValuesTrajFV3JEDI.h"

namespace fv3jedi {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
StateFV3JEDI::StateFV3JEDI(const GeometryFV3JEDI & resol,
                         const oops::Variables & vars,
                         const util::DateTime & vt)
  : fields_(new FieldsFV3JEDI(resol, vars, vt)), stash_()
{
  oops::Log::trace() << "StateFV3JEDI::StateFV3JEDI created." << std::endl;
}
// -----------------------------------------------------------------------------
StateFV3JEDI::StateFV3JEDI(const GeometryFV3JEDI & resol,
                           const eckit::Configuration & file)
  : fields_(), stash_()
{
  const std::vector<std::string> *vv;

  if (file.has("variables"))
    vv = new std::vector<std::string>(file.getStringVector("variables"));
  else
    vv = new std::vector<std::string>({"u", "v", "pt", "delp", "q"});

  oops::Variables vars(*vv);
  fields_.reset(new FieldsFV3JEDI(resol, vars, util::DateTime()));

  if (file.has("analytic_init"))
    fields_->analytic_init(file, resol);
  else
    fields_->read(file);

  ASSERT(fields_);
  oops::Log::trace() << "StateFV3JEDI::StateFV3JEDI created and read in."
                     << std::endl;
}
// -----------------------------------------------------------------------------
StateFV3JEDI::StateFV3JEDI(const GeometryFV3JEDI & resol,
                           const StateFV3JEDI & other)
  : fields_(new FieldsFV3JEDI(*other.fields_, resol)), stash_()
{
  ASSERT(fields_);
  oops::Log::trace() << "StateFV3JEDI::StateFV3JEDI created by interpolation."
                     << std::endl;
}
// -----------------------------------------------------------------------------
StateFV3JEDI::StateFV3JEDI(const StateFV3JEDI & other)
  : fields_(new FieldsFV3JEDI(*other.fields_)), stash_()
{
  ASSERT(fields_);
  oops::Log::trace() << "StateFV3JEDI::StateFV3JEDI copied." << std::endl;
}
// -----------------------------------------------------------------------------
StateFV3JEDI::~StateFV3JEDI() {
  oops::Log::trace() << "StateFV3JEDI::StateFV3JEDI destructed." << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
StateFV3JEDI & StateFV3JEDI::operator=(const StateFV3JEDI & rhs) {
  ASSERT(fields_);
  *fields_ = *rhs.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
/// Get state values at observation locations
// -----------------------------------------------------------------------------
void StateFV3JEDI::getValues(const ioda::Locations & locs,
                             const oops::Variables & vars,
                             ufo::GeoVaLs & cols) const {
  fields_->getValues(locs, vars, cols);
}
// -----------------------------------------------------------------------------
void StateFV3JEDI::getValues(const ioda::Locations & locs,
                             const oops::Variables & vars,
                             ufo::GeoVaLs & cols,
                             const GetValuesTrajFV3JEDI & traj) const {
  fields_->getValues(locs, vars, cols, traj);
}
// -----------------------------------------------------------------------------
/// Interpolate full fields
// -----------------------------------------------------------------------------
void StateFV3JEDI::changeResolution(const StateFV3JEDI & other) {
  fields_->changeResolution(*other.fields_);
  oops::Log::trace() << "StateFV3JEDI changed resolution" << std::endl;
}
// -----------------------------------------------------------------------------
/// Interactions with Increments
// -----------------------------------------------------------------------------
StateFV3JEDI & StateFV3JEDI::operator+=(const IncrementFV3JEDI & dx) {
  ASSERT(this->validTime() == dx.validTime());
  ASSERT(fields_);
  fields_->add(dx.fields());
  return *this;
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void StateFV3JEDI::read(const eckit::Configuration & files) {
  fields_->read(files);
}
// -----------------------------------------------------------------------------
void StateFV3JEDI::analytic_init(const eckit::Configuration & files,
                                 const GeometryFV3JEDI & resol) {
  fields_->analytic_init(files, resol);
}
// -----------------------------------------------------------------------------
void StateFV3JEDI::write(const eckit::Configuration & files) const {
  fields_->write(files);
}
// -----------------------------------------------------------------------------
void StateFV3JEDI::print(std::ostream & os) const {
  os << std::endl << "  Valid time: " << validTime();
  os << *fields_;
}
// -----------------------------------------------------------------------------
/// For accumulator
// -----------------------------------------------------------------------------
void StateFV3JEDI::zero() {
  fields_->zero();
}
// -----------------------------------------------------------------------------
void StateFV3JEDI::accumul(const double & zz, const StateFV3JEDI & xx) {
  fields_->axpy(zz, *xx.fields_);
}
// -----------------------------------------------------------------------------

}  // namespace fv3jedi
