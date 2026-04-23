/*
 * (C) Copyright 2017-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field.h"
#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/ModelData/ModelData.h"

namespace fv3jedi {

class VertRemap : public util::Printable {
 public:
  // Constructor
  VertRemap(const fv3jedi::Geometry & geom, const atlas::FieldSet & fsetOrogNew);

  // Remapping method
  atlas::FieldSet remap(atlas::FieldSet & fset);

  // Constants
  const double beta_    = -6.5e-3;  // lapse rate for temperature extrapolation
  const double epsilon_ = 1.e-9;    // minimum lapse rate so denominators aren't too small

  // Private methods and variables
 private:
  // Methods
  void print(std::ostream &) const;

  // Data
  atlas::Field fieldZsNew_;

  // Model data for Vader
  const ModelData modelData_;
};

}  // namespace fv3jedi
