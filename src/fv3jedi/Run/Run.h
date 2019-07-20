/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3JEDI_RUN_RUN_H_
#define FV3JEDI_RUN_RUN_H_

#include "oops/runs/Run.h"
#include "fv3jedi/Run/Run.interface.h"

namespace fv3jedi {

/*!
 *  Run encapsulates one FV3JEDI/OOPS run.
 */

// -----------------------------------------------------------------------------

class Run : public oops::Run {
 public:
  Run(int, char **);
  ~Run();
};

// -----------------------------------------------------------------------------

}  // namespace fv3jedi
#endif  // FV3JEDI_RUN_RUN_H_
