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

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

#include "oops/mpi/mpi.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "fv3jedi/FieldMetadata/FieldsMetadata.h"
#include "fv3jedi/Geometry/Geometry.interface.h"
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
  static const std::string classname() {return "fv3jedi::Geometry";}

  explicit Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
  Geometry(const Geometry &);
  ~Geometry();

  Geometry & operator=(const Geometry &) = delete;

  bool levelsAreTopDown() const {return true;}

  GeometryIterator begin() const;
  GeometryIterator end() const;
  std::vector<double> verticalCoord(std::string &) const;

  F90geom & toFortran() {return keyGeom_;}
  const F90geom & toFortran() const {return keyGeom_;}
  const eckit::mpi::Comm & getComm() const {return comm_;}
  const atlas::FunctionSpace & functionSpace() const {return functionSpaceIncludingHalo_;}
  atlas::FunctionSpace & functionSpace() {return functionSpaceIncludingHalo_;}
  const atlas::FieldSet & extraFields() const {return extraFields_;}
  atlas::FieldSet & extraFields() {return extraFields_;}
  void latlon(std::vector<double> &, std::vector<double> &, const bool) const;

  std::vector<size_t> variableSizes(const oops::Variables &) const;

  const FieldsMetadata & fieldsMetaData() const {return *fieldsMeta_;}

  // Functions to retrieve geometry features
  const std::vector<double> & ak() const {return ak_;}
  const std::vector<double> & bk() const {return bk_;}
  const double & pTop() const {return pTop_;}
  const int & nLevels() const {return nLevels_;}

 private:
  void print(std::ostream &) const;

  F90geom keyGeom_;
  const eckit::mpi::Comm & comm_;
  atlas::FunctionSpace functionSpace_;
  atlas::FunctionSpace functionSpaceIncludingHalo_;
  atlas::FieldSet extraFields_;
  std::shared_ptr<FieldsMetadata> fieldsMeta_;
  std::vector<double> ak_;
  std::vector<double> bk_;
  int nLevels_;
  double pTop_;
};
// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
