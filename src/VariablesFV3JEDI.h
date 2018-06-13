/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef FV3_JEDI_SRC_VARIABLESFV3JEDI_H_
#define FV3_JEDI_SRC_VARIABLESFV3JEDI_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace fv3jedi {

// -----------------------------------------------------------------------------

class VariablesFV3JEDI : public util::Printable,
                         private util::ObjectCounter<VariablesFV3JEDI> {
 public:
  static const std::string classname() {return "fv3jedi::VariablesFV3JEDI";}

  explicit VariablesFV3JEDI(const oops::Variables &);
  explicit VariablesFV3JEDI(const eckit::Configuration &);

  ~VariablesFV3JEDI();

  VariablesFV3JEDI(const VariablesFV3JEDI &);

  const int * toFortran() const {return &fvars_[0];}

 private:
  void print(std::ostream &) const;
  void setF90(const std::vector<std::string>);
  std::vector<int> fvars_;
};

// -----------------------------------------------------------------------------

}  // namespace fv3jedi

#endif  // FV3_JEDI_SRC_VARIABLESFV3JEDI_H_
