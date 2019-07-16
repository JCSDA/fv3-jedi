/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_GEOMETRY_GEOMETRYFV3JEDI_H_
#define FV3JEDI_GEOMETRY_GEOMETRYFV3JEDI_H_

#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "GeometryFV3JEDIFortran.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {

// -----------------------------------------------------------------------------
/// GeometryFV3JEDI handles geometry for FV3JEDI model.

class GeometryFV3JEDI : public util::Printable,
                      private util::ObjectCounter<GeometryFV3JEDI> {
 public:
  static const std::string classname() {return "fv3jedi::GeometryFV3JEDI";}

  explicit GeometryFV3JEDI(const eckit::Configuration &);
  GeometryFV3JEDI(const GeometryFV3JEDI &);
  ~GeometryFV3JEDI();

  F90geom & toFortran() {return keyGeom_;}
  const F90geom & toFortran() const {return keyGeom_;}

 private:
  GeometryFV3JEDI & operator=(const GeometryFV3JEDI &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // FV3JEDI_GEOMETRY_GEOMETRYFV3JEDI_H_
