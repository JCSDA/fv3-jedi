/*
 * (C) Copyright 2021-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <unordered_set>

#include "fv3jedi/Geometry/Geometry.h"
#include "fv3jedi/IO/Utils/IOBase.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace fv3jedi {

// -------------------------------------------------------------------------------------------------

IOFactory::IOFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    oops::Log::error() << name << " already registered in fv3jedi::IOFactory." << std::endl;
    ABORT("Element already registered in fv3jedi::IOFactory.");
  }
  getMakers()[name] = this;
}

// -------------------------------------------------------------------------------------------------

IOBase * IOFactory::create(const Geometry & geom, const IOParametersBase & params) {
  oops::Log::trace() << "IOBase::create starting" << std::endl;
  const std::string &id = params.filetype.value().value();
  typename std::map<std::string, IOFactory*>::iterator jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    oops::Log::error() << id << " does not exist in fv3jedi::IOFactory." << std::endl;
    ABORT("Element does not exist in fv3jedi::IOFactory.");
  }
  IOBase * ptr = jloc->second->make(geom, params);
  oops::Log::trace() << "IOBase::create done" << std::endl;
  return ptr;
}

// -------------------------------------------------------------------------------------------------

std::unique_ptr<IOParametersBase>
IOFactory::createParameters(const std::string &name) {
  typename std::map<std::string, IOFactory*>::iterator it =
      getMakers().find(name);
  if (it == getMakers().end()) {
    throw std::runtime_error(name + " does not exist in fv3jedi::IOFactory");
  }
  return it->second->makeParameters();
}

// -------------------------------------------------------------------------------------------------

IOBase::IOBase(const Geometry & geom, const eckit::LocalConfiguration conf) {
  oops::Log::trace() << "IOBase::IOBase starting" << std::endl;
  // If conf has 'field io names' then extract from the config and overwrite fieldIoNames_
  if (conf.has("field io names")) {
    fieldIoNames_ = conf.getSubConfiguration("field io names");
  }
  // If conf has 'field io scaling' then extract from the config and overwrite fieldIoScaling_
  if (conf.has("field io scaling")) {
    fieldIoScaling_ = conf.getSubConfiguration("field io scaling");
  }

  // Check the configs for long name correctness
  const std::vector<std::string> fmdLongNames = geom.fieldsMetaData().getLongNames();
  std::unordered_set<std::string> validNames(fmdLongNames.begin(), fmdLongNames.end());

  // Get the keys from the configs
  const std::vector<std::string> fieldIoNamesKeys = fieldIoNames_.keys();
  const std::vector<std::string> fieldIoScalingKeys = fieldIoScaling_.keys();

  // Check for key validity
  for (const std::string& key : fieldIoNamesKeys) {
    const std::string msg = "The \"field io names\" configuration contains \"" + key +
                      "\", which is not part of the field metadata.";
    ASSERT_MSG(validNames.find(key) != validNames.end(), msg);
  }
  for (const std::string& key : fieldIoScalingKeys) {
    const std::string msg = "The \"field io scaling\" configuration contains \"" + key +
                      "\", which is not part of the field metadata.";
    ASSERT_MSG(validNames.find(key) != validNames.end(), msg);
  }

  oops::Log::trace() << "IOBase::IOBase done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void IOBase::readBase(State & x) const {
  // Call read method from the child class
  this->read(x, fieldIoNames_, fieldIoScaling_);
}

// -------------------------------------------------------------------------------------------------

void IOBase::readBase(Increment & dx) const {
  // Call read method from the child class
  this->read(dx, fieldIoNames_, fieldIoScaling_);
}

// -------------------------------------------------------------------------------------------------

void IOBase::writeBase(const State & x) const {
  // Call write method from the child class
  this->write(x, fieldIoNames_, fieldIoScaling_);
}

// -------------------------------------------------------------------------------------------------

void IOBase::writeBase(const Increment & dx) const {
  // Call write method from the child class
  this->write(dx, fieldIoNames_, fieldIoScaling_);
}

// -------------------------------------------------------------------------------------------------

}  // namespace fv3jedi
