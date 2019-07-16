/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Run/RunFV3JEDI.h"
#include "fv3jedi/Utilities/FV3JEDITraits.h"
#include "oops/runs/Forecast.h"

int main(int argc,  char ** argv) {
  fv3jedi::RunFV3JEDI run(argc, argv);
  oops::Forecast<fv3jedi::FV3JEDITraits> fc;
  run.execute(fc);
  return 0;
}
