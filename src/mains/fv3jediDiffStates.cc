/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Run/Run.h"
#include "fv3jedi/Utilities/Traits.h"
#include "oops/runs/DiffStates.h"

int main(int argc,  char ** argv) {
  fv3jedi::Run run(argc, argv);
  oops::DiffStates<fv3jedi::Traits> fc;
  return run.execute(fc);
}
