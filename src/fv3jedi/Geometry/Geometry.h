/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_GEOMETRY_GEOMETRY_H_
#define FV3JEDI_GEOMETRY_GEOMETRY_H_

#include <ostream>
#include <string>

#include "eckit/mpi/Comm.h"

#include "oops/parallel/mpi/mpi.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/Geometry/Geometry.interface.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {

// -----------------------------------------------------------------------------
/// Geometry handles geometry for FV3JEDI model.

class Geometry : public util::Printable,
                      private util::ObjectCounter<Geometry> {
 public:
  static const std::string classname() {return "fv3jedi::Geometry";}

  explicit Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
  Geometry(const Geometry &);
  ~Geometry();

  F90geom & toFortran() {return keyGeom_;}
  const F90geom & toFortran() const {return keyGeom_;}
  const eckit::mpi::Comm & getComm() const {return comm_;}

 private:
  Geometry & operator=(const Geometry &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
  const eckit::mpi::Comm & comm_;
};
// -----------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // FV3JEDI_GEOMETRY_GEOMETRY_H_
