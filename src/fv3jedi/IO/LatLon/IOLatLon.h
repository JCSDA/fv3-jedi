/*
 * (C) Copyright 2021 UCAR.
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
#include "IOLatLon.interface.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

class IOLatLonParameters : public IOParametersBase {
  OOPS_CONCRETE_PARAMETERS(IOLatLonParameters, IOParametersBase)

 public:
  // Path and name of the file to be written
  oops::Parameter<std::string> filename{"filename", "filename to be written",
                                        "Data/fv3jedi.latlon.", this};

  // Optionally config may contain member
  oops::OptionalParameter<int> member{"member", "ensemble member number", this};

  // Things BUMP puts in write config not needed in IO but specified to avoid failures
  oops::OptionalParameter<std::string> bumpparameter{"parameter", "bump parameter", this};
  oops::OptionalParameter<util::DateTime> bumpdatetime{"date", "bump datetime", this};
};

// -------------------------------------------------------------------------------------------------
class IOLatLon : public IOBase, private util::ObjectCounter<IOLatLon> {
 public:
  static const std::string classname() {return "fv3jedi::IOLatLon";}

  typedef IOLatLonParameters Parameters_;

  IOLatLon(const Geometry &, const Parameters_ &);
  ~IOLatLon();
  void read(State &) const override;
  void read(Increment &) const override;
  void write(const State &) const override;
  void write(const Increment &) const override;

 private:
  F90iolatlon objectKeyForFortran_;
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
