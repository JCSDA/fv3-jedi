/*
 * (C) Copyright 2018 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "FV3JEDITraits.h"
#include "instantiateFV3JEDIVarChangeFactory.h"
#include "oops/runs/EstimateParams.h"
#include "RunFV3JEDI.h"

int main(int argc,  char ** argv) {
  fv3jedi::RunFV3JEDI run(argc, argv);
  fv3jedi::instantiateFV3JEDIVarChangeFactory();
  oops::EstimateParams<fv3jedi::FV3JEDITraits> dir;
  run.execute(dir);
  return 0;
}
