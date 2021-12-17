/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#pragma once

#include <algorithm>
#include <iterator>
#include <map>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/config/YAMLConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

#include "fv3jedi/FieldMetadata/FieldsMetadataParameters.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {

  // -----------------------------------------------------------------------------------------------

  class FieldMetadata : public util::Printable {
   public:
    explicit FieldMetadata(const std::string fieldIOName) {
      fieldIOName_ = fieldIOName;
    }

    // FieldIOName
    std::string getFieldIOName() const {return fieldIOName_;}
    void setFieldIOName(std::string fieldIOName) {fieldIOName_ = fieldIOName;}

    // FieldName
    std::string getFieldName() const {return fieldName_;}
    void setFieldName(std::string fieldName) {fieldName_ = fieldName;}

    // Kind
    std::string getKind() const {return kind_;}
    void setKind(std::string kind) {kind_ = kind;}

    // Levels
    int getLevels() const {return levels_;}
    void setLevels(int levels) {levels_ = levels;}

    // LongName
    std::string getLongName() const {return longName_;}
    void setLongName(std::string longName) {longName_ = longName;}

    // Space
    std::string getSpace() const {return space_;}
    void setSpace(std::string space) {space_ = space;}

    // StaggerLoc
    std::string getStaggerLoc() const {return staggerLoc_;}
    void setStaggerLoc(std::string staggerLoc) {staggerLoc_ = staggerLoc;}

    // Tracer
    bool getTracer() const {return tracer_;}
    void setTracer(bool tracer) {tracer_ = tracer;}

    // Units
    std::string getUnits() const {return units_;}
    void setUnits(std::string units) {units_ = units;}

    // InterpType
    std::string getInterpType() const {return interpType_;}
    void setInterpType(std::string interpType) {interpType_ = interpType;}

    // IO file
    std::string getIOFile() const {return io_file_;}
    void setIOFile(std::string io_file) {io_file_ = io_file;}

    // Validity check on kind
    void checkKindValid(const std::string fieldIOName, const std::string kind) const {
      auto result = std::find(kindVal_.begin(), kindVal_.end(), kind);
      if (result == std::end(kindVal_)) {
        oops::Log::debug() << "Key: " << fieldIOName << " failed due to invalid kind: "
                           << kind << ". Options include: " << kindVal_ << std::endl;
        ABORT("FieldMetadata::checkKindValid failed, run again with debug prints");
      }
    }

    // Validity check on stagger location
    void checkStaggerLocValid(const std::string fieldIOName, const std::string staggerLoc) const {
      auto result = std::find(staggerLocVal_.begin(), staggerLocVal_.end(), staggerLoc);
      if (result == std::end(staggerLocVal_)) {
        oops::Log::debug() << "Key: " << fieldIOName << " failed due to invalid stagger location: "
                           << staggerLoc << ". Options include: " << staggerLocVal_ << std::endl;
        ABORT("FieldMetadata::checkStaggerLocValid failed, run again with debug prints");
      }
    }

    // Validity check on space
    void checkSpaceValid(const std::string fieldIOName, const std::string space) const {
      auto result = std::find(spaceVal_.begin(), spaceVal_.end(), space);
      if (result == std::end(spaceVal_)) {
        oops::Log::debug() << "Key: " << fieldIOName << " failed due to invalid space: "
                           << space << ". Options include: " << spaceVal_ << std::endl;
        ABORT("FieldMetadata::checkSpace failed, run again with debug prints");
      }
    }

    // Validity check on interpolation type
    void checkInterpTypeValid(const std::string fieldIOName, const std::string interpType) const {
      auto result = std::find(interpTypeVal_.begin(), interpTypeVal_.end(), interpType);
      if (result == std::end(interpTypeVal_)) {
        oops::Log::debug() << "Key: " << fieldIOName << " failed due to invalid interType: "
                           << interpType << ". Options include: " << interpTypeVal_ << std::endl;
        ABORT("FieldMetadata::checkSpace failed, run again with debug prints");
      }
    }

    // Validity check on string level
    void checkLevelValid(const std::string fieldIOName, const std::string level) const {
      auto result = std::find(levelVal_.begin(), levelVal_.end(), level);
      if (result == std::end(levelVal_) && is_number(level) == false) {
        oops::Log::debug() << "Key: " << fieldIOName << " failed due to invalid level choice: "
                           << level << ". Options include: " << levelVal_ << "or an integer"
                           << std::endl;
        ABORT("FieldMetadata::checkLevelValid failed, run again with debug prints");
      }
    }

    // Function to check if string is a number
    bool is_number(const std::string& s) const {
      return !s.empty() && std::find_if(s.begin(), s.end(), [](unsigned char c)
                                        { return !std::isdigit(c); }) == s.end();
    }

   private:
    std::string fieldIOName_;
    std::string fieldName_;
    std::string kind_;
    int levels_;
    std::string longName_;
    std::string staggerLoc_;
    std::string space_;
    bool tracer_;
    std::string units_;
    std::string interpType_;
    std::string io_file_;

    const std::vector<std::string> kindVal_ = {"double", "integer"};
    const std::vector<std::string> levelVal_ = {"full", "half"};
    const std::vector<std::string> spaceVal_ = {"vector", "magnitude", "direction"};
    const std::vector<std::string> staggerLocVal_ = {"center", "eastwest", "northsouth", "corner"};
    const std::vector<std::string> interpTypeVal_ = {"integer", "nearest", "default"};

    void print(std::ostream & os) const {
      os << " FieldIOName: " << fieldIOName_ << "\n";
      os << " FieldName: " << fieldName_ << "\n";
      os << " Kind: " << kind_ << "\n";
      os << " Levels: " << levels_ << "\n";
      os << " LongName: " << longName_ << "\n";
      os << " Space: " << space_ << "\n";
      os << " StaggerLocation: " << staggerLoc_ << "\n";
      os << " Tracer: " << tracer_ << "\n";
      os << " Units: " << units_ << "\n";
      os << " InterpType: " << interpType_ << "\n";
      os << " IOFile: " << io_file_ << "\n";
    }
  };

  // -----------------------------------------------------------------------------------------------

  class FieldsMetadata {
   public:
    typedef FieldsMetadataParameters Parameters_;
    FieldsMetadata(const Parameters_ &, int &);

    FieldMetadata getField(const std::string &) const;

    // Function to return number of levels given the longName
    size_t getLevelsFromLongName(const std::string &) const;

    // Get long name from IO name
    oops::Variables LongNameFromIONameLongNameOrFieldName(const oops::Variables &) const;

   private:
    std::map<std::string, FieldMetadata> fields_;
  };

  // -----------------------------------------------------------------------------------------------

}  // namespace fv3jedi
