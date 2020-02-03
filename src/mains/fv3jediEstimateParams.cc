/*
 * (C) Copyright 2018 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Run/Run.h"
#include "fv3jedi/Utilities/Traits.h"
#include "fv3jedi/VariableChanges/instantiateVarChangeFactories.h"
#include "oops/runs/EstimateParams.h"

int main(int argc,  char ** argv) {
  fv3jedi::Run run(argc, argv);
  fv3jedi::instantiateVarChangeFactories();
  oops::EstimateParams<fv3jedi::Traits> dir;
  run.execute(dir);
  return 0;
}
