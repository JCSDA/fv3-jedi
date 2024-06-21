/*
 * (C) Copyright 2023 UCAR.
 * (C) Crown copyright 2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Utilities/Traits.h"
#include "oops/generic/instantiateLinearModelFactory.h"
#include "oops/generic/instantiateModelFactory.h"
#include "oops/runs/LinearizationError.h"
#include "oops/runs/Run.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  oops::instantiateLinearModelFactory<fv3jedi::Traits>();
  oops::instantiateModelFactory<fv3jedi::Traits>();
  oops::LinearizationError<fv3jedi::Traits> linearizationError;
  return run.execute(linearizationError);
}
