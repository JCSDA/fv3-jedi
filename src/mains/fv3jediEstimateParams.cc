/*
 * (C) Copyright 2018 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Run/RunFV3JEDI.h"
#include "fv3jedi/Utilities/FV3JEDITraits.h"
#include "fv3jedi/Utilities/instantiateFV3JEDIVarChangeFactories.h"
#include "oops/runs/EstimateParams.h"

int main(int argc,  char ** argv) {
  fv3jedi::RunFV3JEDI run(argc, argv);
  fv3jedi::instantiateFV3JEDIVarChangeFactories();
  oops::EstimateParams<fv3jedi::FV3JEDITraits> dir;
  run.execute(dir);
  return 0;
}
