/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "fv3jedi/Utilities/Traits.h"

#include "oops/runs/Run.h"
#include "test/interface/ObsLocalization.h"
#include "fv3jedi/ObsLocalization/instantiateObsLocFactory.h"
#include "ufo/ObsTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  fv3jedi::instantiateObsLocFactory();
  test::ObsLocalization<fv3jedi::Traits, ufo::ObsTraits> tests;
  return run.execute(tests);
}
