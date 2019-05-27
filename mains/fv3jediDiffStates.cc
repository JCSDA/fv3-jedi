/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "FV3JEDITraits.h"
#include "oops/runs/DiffStates.h"
#include "RunFV3JEDI.h"

int main(int argc,  char ** argv) {
  fv3jedi::RunFV3JEDI run(argc, argv);
  oops::DiffStates<fv3jedi::FV3JEDITraits> fc;
  run.execute(fc);
  return 0;
}
