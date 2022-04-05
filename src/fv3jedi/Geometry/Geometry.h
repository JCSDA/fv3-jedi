/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

#include "oops/mpi/mpi.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/FieldMetadata/FieldsMetadata.h"
#include "fv3jedi/Geometry/Geometry.interface.h"
#include "fv3jedi/Geometry/GeometryParameters.h"
#include "fv3jedi/GeometryIterator/GeometryIterator.h"

namespace oops {
  class Variables;
}

namespace fv3jedi {
  class GeometryIterator;

// -------------------------------------------------------------------------------------------------
/// Geometry handles geometry for FV3JEDI model.

class Geometry : public util::Printable,
                 private util::ObjectCounter<Geometry> {
 public:
  typedef GeometryParameters Parameters_;
  static const std::string classname() {return "fv3jedi::Geometry";}

  explicit Geometry(const Parameters_ &, const eckit::mpi::Comm &);
  Geometry(const Geometry &);
  ~Geometry();

  GeometryIterator begin() const;
  GeometryIterator end() const;
  std::vector<double> verticalCoord(std::string &) const;

  F90geom & toFortran() {return keyGeom_;}
  const F90geom & toFortran() const {return keyGeom_;}
  const eckit::mpi::Comm & getComm() const {return comm_;}
  atlas::FunctionSpace * atlasFunctionSpace() const {return atlasFunctionSpace_.get();}
  atlas::FieldSet * atlasFieldSet() const {return atlasFieldSet_.get();}
  void latlon(std::vector<double> &, std::vector<double> &, const bool) const;

  std::vector<size_t> variableSizes(const oops::Variables &) const;

  const FieldsMetadata & fieldsMetaData() const {return *fieldsMeta_;}

 private:
  Geometry & operator=(const Geometry &);
  void print(std::ostream &) const;
  F90geom keyGeom_;
  const eckit::mpi::Comm & comm_;
  std::unique_ptr<atlas::functionspace::PointCloud> atlasFunctionSpace_;
  std::unique_ptr<atlas::functionspace::PointCloud> atlasFunctionSpaceIncludingHalo_;
  std::unique_ptr<atlas::FieldSet> atlasFieldSet_;
  std::shared_ptr<FieldsMetadata> fieldsMeta_;
};
// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
