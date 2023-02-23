/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Geometry/TimeInvariantFieldsHelpers.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"

namespace fv3jedi {

extern "C" {
  void fv3jedi_nominal_surface_pressure_f90(atlas::field::FieldSetImpl *);
}

void insertDerivedTimeInvariantFields(atlas::FieldSet & fset,
                                      const oops::Variables & varsToAdd) {
  for (const auto & var : varsToAdd.variables()) {
    if (var == "nominal_surface_pressure") {
      // Check requested field doesn't already exist
      ASSERT(!fset.has(var));
      // Check necessary inputs are present
      ASSERT(fset.has("filtered_orography"));

      // Add new atlas::Field for nominal surface pressure
      const auto & orog = fset.field("filtered_orography");
      atlas::Field nsp =
          orog.functionspace().template createField<double>(atlas::option::name(var)
                                                            | atlas::option::levels(1));
      fset.add(nsp);

      // Call fortran to set values in the new NSP field
      fv3jedi_nominal_surface_pressure_f90(fset.get());
    } else {
      throw eckit::UserError("Unsupported time-invariant derived field: " + var);
    }
  }
}

}  // namespace fv3jedi
