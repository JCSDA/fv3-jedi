/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_GEOMETRY_GEOMETRY_H_
#define FV3JEDI_GEOMETRY_GEOMETRY_H_

#include <memory>
#include <ostream>
#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

#include "oops/parallel/mpi/mpi.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/FieldMetadata/FieldsMetadata.h"
#include "fv3jedi/Geometry/Geometry.interface.h"
#include "fv3jedi/GeometryIterator/GeometryIterator.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {

  class GeometryIterator;

// -------------------------------------------------------------------------------------------------
/// Geometry handles geometry for FV3JEDI model.

class Geometry : public util::Printable,
                 private util::ObjectCounter<Geometry> {
 public:
  static const std::string classname() {return "fv3jedi::Geometry";}

  explicit Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
  Geometry(const Geometry &);
  ~Geometry();

  GeometryIterator begin() const;
  GeometryIterator end() const;

  F90geom & toFortran() {return keyGeom_;}
  const F90geom & toFortran() const {return keyGeom_;}
  const eckit::mpi::Comm & getComm() const {return comm_;}
  atlas::FunctionSpace * atlasFunctionSpace() const {return atlasFunctionSpace_.get();}
  atlas::FieldSet * atlasFieldSet() const {return atlasFieldSet_.get();}

 private:
  Geometry & operator=(const Geometry &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
  const eckit::mpi::Comm & comm_;
  std::unique_ptr<atlas::functionspace::PointCloud> atlasFunctionSpace_;
  std::unique_ptr<atlas::FieldSet> atlasFieldSet_;
  const FieldsMetadata fieldsMeta_;
};
// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // FV3JEDI_GEOMETRY_GEOMETRY_H_
