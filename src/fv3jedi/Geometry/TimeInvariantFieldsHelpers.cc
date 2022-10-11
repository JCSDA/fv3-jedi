/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Geometry/TimeInvariantFieldsHelpers.h"

#include "atlas/field.h"
#include "atlas/functionspace/FunctionSpace.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"

namespace detail {

void computeExample(atlas::Field & exampleOutput,
                    const atlas::Field & exampleInput0,
                    const atlas::Field & exampleInput1) {
  // math here
}

}  // namespace detail

namespace fv3jedi {

void insertDerivedTimeInvariantFields(atlas::FieldSet & fset, const oops::Variables & varsToAdd) {
  for (const auto & var : varsToAdd.variables()) {
    if (var == "example_output") {
      // Check requested field doesn't already exist
      ASSERT(!fset.has(var));
      // Check necessary inputs are present
      ASSERT(fset.has("example_input0"));
      ASSERT(fset.has("example_input1"));

      const auto & exampleInput0 = fset.field("example_input0");
      const auto & exampleInput1 = fset.field("example_input1");
      atlas::Field exampleOutput =
          exampleInput0.functionspace().template createField<double>(atlas::option::name(var)
                                                                     | atlas::option::levels(1));
      detail::computeExample(exampleOutput, exampleInput0, exampleInput1);
      fset.add(exampleOutput);
    } else {
      throw eckit::UserError("Unsupported time-invariant derived field: " + var);
    }
  }
}

}  // namespace fv3jedi
