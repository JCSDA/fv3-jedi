/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Utilities/Traits.h"
#include "ioda/instantiateObsLocFactory.h"
#include "oops/runs/LETKF.h"
#include "oops/runs/Run.h"
#include "ufo/instantiateObsFilterFactory.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ioda::instantiateObsLocFactory<fv3jedi::Traits>();
  ufo::instantiateObsFilterFactory<fv3jedi::Traits>();
  oops::LETKF<fv3jedi::Traits> letkf;
  run.execute(letkf);
  return 0;
}
