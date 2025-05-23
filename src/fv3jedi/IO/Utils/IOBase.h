/*
 * (C) Copyright 2021-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/Variables.h"
#include "oops/util/AssociativeContainers.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/OptionalPolymorphicParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/Printable.h"

namespace fv3jedi {
  class Geometry;
  class Increment;
  class State;

// -------------------------------------------------------------------------------------------------

class IOBase : public util::Printable, private boost::noncopyable {
 public:
  explicit IOBase(const Geometry &, const eckit::LocalConfiguration);
  virtual ~IOBase() {}

  void readBase(State &) const;
  void readBase(Increment &) const;
  void writeBase(const State &) const;
  void writeBase(const Increment &) const;

 private:
  // Child read/write methods
  virtual void read(State &, const eckit::LocalConfiguration &,
                    const eckit::LocalConfiguration &) const = 0;
  virtual void read(Increment &, const eckit::LocalConfiguration &,
                    const eckit::LocalConfiguration &) const = 0;
  virtual void write(const State &, const eckit::LocalConfiguration &,
                     const eckit::LocalConfiguration &) const = 0;
  virtual void write(const Increment &, const eckit::LocalConfiguration &,
                     const eckit::LocalConfiguration &) const = 0;

  // Child print method
  virtual void print(std::ostream &) const = 0;

  // Configuration holding the field names as used in the files
  eckit::LocalConfiguration fieldIoNames_;
  // Configuration holding scaling factors used to transition units between forecast model
  // and JEDI
  eckit::LocalConfiguration fieldIoScaling_;
};

// -------------------------------------------------------------------------------------------------

class IOParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(IOParametersBase, Parameters)
 public:
  oops::OptionalParameter<std::string> type{"type", this};
  oops::OptionalParameter<std::string> exp{"exp", this};
  oops::OptionalParameter<int> member{"member", this};
  oops::OptionalParameter<std::string> memberPattern{"member pattern", this};
  oops::OptionalParameter<util::DateTime> date{"date", this};
  oops::OptionalParameter<int> iteration{"iteration", this};
  oops::OptionalParameter<std::string> prefix{"prefix", this};
  oops::Parameter<bool> dateCols{"date colons", true, this};
  oops::OptionalParameter<std::string> filetype{"filetype", this};
  oops::OptionalParameter<eckit::LocalConfiguration> fieldIoNames{"field io names", this};
};

// -------------------------------------------------------------------------------------------------

class IOFactory;

// -------------------------------------------------------------------------------------------------

class IOParametersWrapper : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(IOParametersWrapper, Parameters)
 public:
  oops::OptionalPolymorphicParameter<IOParametersBase, IOFactory> ioParameters{"filetype", this};
};

// -------------------------------------------------------------------------------------------------

class IOFactory {
 public:
  static IOBase * create(const Geometry &, const IOParametersBase &params);

  static std::unique_ptr<IOParametersBase> createParameters(const std::string &name);

  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~IOFactory() = default;

 protected:
  explicit IOFactory(const std::string &name);

 private:
  virtual IOBase * make(const Geometry &, const IOParametersBase &) = 0;

  virtual std::unique_ptr<IOParametersBase> makeParameters() const = 0;

  static std::map < std::string, IOFactory * > & getMakers() {
    static std::map < std::string, IOFactory * > makers_;
    return makers_;
  }
};

// -------------------------------------------------------------------------------------------------

template<class T>
class IOMaker : public IOFactory {
  typedef typename T::Parameters_ Parameters_;

  IOBase * make(const Geometry & geom, const IOParametersBase & params) override {
    return new T(geom, dynamic_cast<const Parameters_&>(params));
  }

  std::unique_ptr<IOParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit IOMaker(const std::string & name) : IOFactory(name) {}
};

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
