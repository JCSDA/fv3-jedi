/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "FV3JEDITraits.h"
#include "instantiateLocalizationFactory.h"
#include "instantiateFV3JEDIVarChangeFactory.h"
#include "oops/runs/Variational.h"
#include "RunFV3JEDI.h"

int main(int argc,  char ** argv) {
  fv3jedi::RunFV3JEDI run(argc, argv);
  fv3jedi::instantiateLocalizationFactory();
  fv3jedi::instantiateFV3JEDIVarChangeFactory();
  oops::Variational<fv3jedi::FV3JEDITraits> var;
  run.execute(var);
  return 0;
}
