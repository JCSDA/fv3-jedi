/*
 * (C) Copyright 2017-2021 UCAR
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

#include "oops/base/LocalIncrement.h"
#include "oops/base/WriteParametersBase.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "fv3jedi/Increment/Increment.interface.h"
#include "fv3jedi/IO/Utils/IOBase.h"

namespace oops {
  class Variables;
}

namespace fv3jedi {
  class Geometry;
  class ModelBiasIncrement;
  class State;

// -------------------------------------------------------------------------------------------------

class DiracParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DiracParameters, Parameters)
 public:
  // Number of Diracs
  oops::RequiredParameter<int> ndir{"ndir", this};
  // Vector (length n) with index in x direction
  oops::RequiredParameter<std::vector<int>> ixdir{"ixdir", this};
  // Vector (length n) with index in y direction
  oops::RequiredParameter<std::vector<int>> iydir{"iydir", this};
  // Vector (length n) with index of level
  oops::RequiredParameter<std::vector<int>> ildir{"ildir", this};
  // Vector (length n) with index of tile
  oops::RequiredParameter<std::vector<int>> itdir{"itdir", this};
  // Vector (length n) with name of field
  oops::RequiredParameter<std::vector<std::string>> ifdir{"ifdir", this};
};

// -------------------------------------------------------------------------------------------------

class IncrementReadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(IncrementReadParameters, Parameters)
 public:
  IOParametersWrapper ioParametersWrapper{this};
  oops::OptionalParameter<bool> setdatetime{"set datetime on read", this};
  oops::OptionalParameter<util::DateTime> datetime{"datetime", this};
};

// -------------------------------------------------------------------------------------------------

class IncrementWriteParameters : public oops::WriteParametersBase {
  OOPS_CONCRETE_PARAMETERS(IncrementWriteParameters, WriteParametersBase)
 public:
  IOParametersWrapper ioParametersWrapper{this};
};

// -------------------------------------------------------------------------------------------------

class Increment : public util::Printable,
                  private util::ObjectCounter<Increment> {
 public:
  static const std::string classname() {return "fv3jedi::Increment";}

  typedef DiracParameters          DiracParameters_;
  typedef IncrementReadParameters  ReadParameters_;
  typedef IncrementWriteParameters WriteParameters_;

/// Constructor, destructor
  Increment(const Geometry &, const oops::Variables &, const util::DateTime &);
  Increment(const Geometry &, const Increment &);
  Increment(const Increment &, const bool);
  virtual ~Increment();

/// Basic operators
  void diff(const State &, const State &);
  void zero();
  void zero(const util::DateTime &);
  void ones();
  Increment & operator =(const Increment &);
  Increment & operator+=(const Increment &);
  Increment & operator-=(const Increment &);
  Increment & operator*=(const double &);
  void axpy(const double &, const Increment &, const bool check = true);
  double dot_product_with(const Increment &) const;
  void schur_product_with(const Increment &);
  void random();
  void dirac(const DiracParameters_ &);


/// Get/Set increment values at grid points
  oops::LocalIncrement getLocal(const GeometryIterator &) const;
  void setLocal(const oops::LocalIncrement &, const GeometryIterator &);

/// ATLAS
  void setAtlas(atlas::FieldSet *) const;
  void toAtlas(atlas::FieldSet *) const;
  void fromAtlas(atlas::FieldSet *);

/// I/O and diagnostics
  void read(const ReadParameters_ &);
  void write(const WriteParameters_ &) const;
  double norm() const;

// Add or remove fields
  void updateFields(const oops::Variables &);

  void updateTime(const util::Duration & dt) {time_ += dt;}

/// Other
  void accumul(const double &, const State &);

/// Serialize and deserialize
  size_t serialSize() const;
  void serialize(std::vector<double> &) const;
  void deserialize(const std::vector<double> &, size_t &);

// Utilities
  std::shared_ptr<const Geometry> geometry() const {return geom_;}
  const oops::Variables & variables() const {return vars_;}
  const oops::Variables & variablesLongName() const {return varsLongName_;}

  const util::DateTime & time() const {return time_;}
  util::DateTime & time() {return time_;}
  const util::DateTime & validTime() const {return time_;}
  util::DateTime & validTime() {return time_;}

  int & toFortran() {return keyInc_;}
  const int & toFortran() const {return keyInc_;}

// Private methods and variables
 private:
  void print(std::ostream &) const;
  F90inc keyInc_;
  std::shared_ptr<const Geometry> geom_;
  oops::Variables vars_;
  oops::Variables varsLongName_;
  util::DateTime time_;
};
// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
