#ifndef FV3JEDI_RUNFV3JEDI_H_
#define FV3JEDI_RUNFV3JEDI_H_

#include "oops/runs/Run.h"

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
#endif  // FV3JEDI_RUNFV3JEDI_H_
