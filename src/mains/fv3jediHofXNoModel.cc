/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Run/Run.h"
#include "fv3jedi/Utilities/Traits.h"
#include "oops/runs/HofXNoModel.h"
#include "ufo/instantiateObsFilterFactory.h"

int main(int argc,  char ** argv) {
  fv3jedi::Run run(argc, argv);
  ufo::instantiateObsFilterFactory<fv3jedi::Traits>();
  oops::HofXNoModel<fv3jedi::Traits> hofx;
  return run.execute(hofx);
}
