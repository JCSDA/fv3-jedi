/*
 * (C) Copyright 2022 NOAA
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "oops/util/DateTime.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "fv3jedi/IO/Utils/IOBase.h"
#include "IOAuxGrid.interface.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

class IOAuxGridParameters : public IOParametersBase {
  OOPS_CONCRETE_PARAMETERS(IOAuxGridParameters, IOParametersBase)

 public:
  // Type of gridded output
  oops::Parameter<std::string> gridtype{"gridtype", "gridtype to be populated",
                                        "gaussian", this};
  // Path and name of the file to be written
  oops::Parameter<std::string> filename{"filename", "filename to be written",
                                        "Data/fv3jedi.auxgrid.", this};

  // Optionally config may contain member
  oops::OptionalParameter<int> member{"member", "ensemble member number", this};
};

// -------------------------------------------------------------------------------------------------
class IOAuxGrid : public IOBase, private util::ObjectCounter<IOAuxGrid> {
 public:
  static const std::string classname() {return "fv3jedi::IOAuxGrid";}

  typedef IOAuxGridParameters Parameters_;

  IOAuxGrid(const Geometry &, const Parameters_ &);
  ~IOAuxGrid();
  void read(State &) const override;
  void read(Increment &) const override;
  void write(const State &) const override;
  void write(const Increment &) const override;

 private:
  F90ioauxgrid objectKeyForFortran_;
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
