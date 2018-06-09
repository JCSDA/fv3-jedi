/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "FV3JEDITraits.h"
#include "RunFV3JEDI.h"
#include "test/interface/Localization.h"

int main(int argc,  char ** argv) {
  fv3jedi::RunFV3JEDI run(argc, argv);
  test::Localization<fv3jedi::FV3JEDITraits> tests;
  run.execute(tests);
  return 0;
}

