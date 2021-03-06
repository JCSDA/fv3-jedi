/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_GEOMETRYITERATOR_GEOMETRYITERATOR_H_
#define FV3JEDI_GEOMETRYITERATOR_GEOMETRYITERATOR_H_

#include <iterator>
#include <string>

#include "eckit/geometry/Point2.h"

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/GeometryIterator/GeometryIterator.interface.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace fv3jedi {

class Geometry;

// -----------------------------------------------------------------------------
class GeometryIterator: public std::iterator<std::forward_iterator_tag,
                                               eckit::geometry::Point2>,
                          public util::Printable,
                          private util::ObjectCounter<GeometryIterator> {
 public:
  static const std::string classname() {return "fv3jedi::GeometryIterator";}

  GeometryIterator(const GeometryIterator &);
  explicit GeometryIterator(const Geometry & geom,
                            const int & iindex = 1, const int & jindex = 1);
  ~GeometryIterator();

  bool operator==(const GeometryIterator &) const;
  bool operator!=(const GeometryIterator &) const;
  eckit::geometry::Point2 operator*() const;
  GeometryIterator& operator++();

// Utilities
  F90iter & toFortran() {return keyIter_;}
  const F90iter & toFortran() const {return keyIter_;}

 private:
  void print(std::ostream &) const;
  F90iter keyIter_;
  const  Geometry & geom_;
};

}  // namespace fv3jedi

#endif  // FV3JEDI_GEOMETRYITERATOR_GEOMETRYITERATOR_H_
