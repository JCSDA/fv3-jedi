/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include <boost/ptr_container/ptr_vector.hpp>

#include "eckit/config/Configuration.h"

#include "oops/util/Printable.h"

#include "vader/vader.h"

#include "fv3jedi/FieldMetadata/FieldsMetadata.h"
#include "fv3jedi/VariableChange/Base/VariableChangeBase.h"

namespace fv3jedi {
  class Geometry;
  class State;

// -------------------------------------------------------------------------------------------------

class VariableChange : public util::Printable {
 public:
  static const std::string classname() {return "fv3jedi::VariableChange";}

  explicit VariableChange(const eckit::Configuration &, const Geometry &);
  ~VariableChange();

  // The bool `force_varchange` is used in the "model-to-analysis" var change applied in the model
  // finalize call; this forces a D-to-A conversion of the winds, even if the model and the
  // analysis both contain A-grid winds and the var change looks like it can be skipped
  void changeVar(State &, const oops::Variables &) const;
  void changeVarInverse(State &, const oops::Variables &, bool force_varchange = false) const;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<VariableChangeBase> variableChange_;
  FieldsMetadata fieldsMetadata_;
  bool run_fv3jedi_;
  bool run_vader_;
  std::unique_ptr<vader::Vader> vader_;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
