/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_RUN_RUNFV3JEDI_H_
#define FV3JEDI_RUN_RUNFV3JEDI_H_

#include "oops/runs/Run.h"
#include "RunFV3JEDIFortran.h"

namespace fv3jedi {

/*!
 *  RunFV3JEDI encapsulates one FV3JEDI/OOPS run.
 */

// -----------------------------------------------------------------------------

class RunFV3JEDI : public oops::Run {
 public:
  RunFV3JEDI(int, char **);
  ~RunFV3JEDI();
};

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_RUN_RUNFV3JEDI_H_
