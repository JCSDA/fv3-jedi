/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/config/Configuration.h"
#include "fv3jedi/GeometryIterator/GeometryIterator.h"
#include "fv3jedi/GeometryIterator/GeometryIterator.interface.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------

namespace fv3jedi {


// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const GeometryIterator& iter):
  geom_(iter.geom_)
{
  fv3jedi_geom_iter_clone_f90(keyIter_, iter.toFortran());
}

// -----------------------------------------------------------------------------

GeometryIterator::GeometryIterator(const Geometry & geom,
                                       const int & iindex, const int & jindex):
  geom_(geom)
{
  fv3jedi_geom_iter_setup_f90(keyIter_, geom_.toFortran(), iindex, jindex);
}


// -----------------------------------------------------------------------------

GeometryIterator::~GeometryIterator() {
  fv3jedi_geom_iter_delete_f90(keyIter_);
}

// -----------------------------------------------------------------------------

bool GeometryIterator::operator==(const GeometryIterator & other) const {
  int equals = 0;
  fv3jedi_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return (equals == 1);
}

// -----------------------------------------------------------------------------

bool GeometryIterator::operator!=(const GeometryIterator & other) const {
  int equals = 0;
  fv3jedi_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return (equals == 0);
}

// -----------------------------------------------------------------------------

eckit::geometry::Point2 GeometryIterator::operator*() const {
  double lat, lon;
  fv3jedi_geom_iter_current_f90(keyIter_, lon, lat);
  return eckit::geometry::Point2(lon, lat);
}

// -----------------------------------------------------------------------------

GeometryIterator& GeometryIterator::operator++() {
  fv3jedi_geom_iter_next_f90(keyIter_);
  return *this;
}

// -----------------------------------------------------------------------------

double GeometryIterator::getOrography() const {
  double orography;
  fv3jedi_geom_iter_orography_f90(keyIter_, orography);
  return orography;
}

// -----------------------------------------------------------------------------

void GeometryIterator::print(std::ostream & os) const {
  double lat, lon;
  fv3jedi_geom_iter_current_f90(keyIter_, lon, lat);
  os << "GeometryIterator, lat/lon: " << lat << " / " << lon << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
