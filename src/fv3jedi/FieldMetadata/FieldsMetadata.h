/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
*/

#pragma once

#include <algorithm>
#include <iostream>
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
#include "oops/util/Printable.h"

#include "fv3jedi/FieldMetadata/FieldsMetadataParameters.h"

namespace eckit {
  class Configuration;
}

namespace fv3jedi {

  // -----------------------------------------------------------------------------------------------

  class FieldMetadata : public util::Printable {
   public:
    explicit FieldMetadata(const std::string longName, int nlev) {
      longName_ = longName;
      nlev_ = nlev;
    }

    // Get functions
    // -------------
    bool getIsTracer() const {return isTracer_;}
    int getNumLevls() const {return numLevls_;}
    std::string getLongName() const {return longName_;}
    std::string getShrtName() const {return shrtName_;}
    std::string getDataKind() const {return dataKind_;}
    std::string getStagrLoc() const {return stagrLoc_;}
    std::string getMathSpac() const {return mathSpac_;}
    std::string getVarUnits() const {return varUnits_;}
    std::string getInOuName() const {return inOuName_;}
    std::string getInOuFile() const {return inOuFile_;}
    std::string getIntrpTyp() const {return intrpTyp_;}
    std::string getIntrpMsk() const {return intrpMsk_;}

    // Whether a field is a specific to the FV3 dycore or to a model using the FV3 dycore; we lump
    // these together as "interface specific" fields. This is contrast to fields that have meaning
    // outside the context of the FV3 dycore or a model, which are used outside of this model
    // interface by JEDI.
    //
    // For now, we determine if a field is interface specific algorithmically:
    // => a field is interface specific IFF (is it staggered OR its long name ends in _cold)
    //
    // But if more fine-grained control is needed, this could be set by a data member, as is done
    // for other metadata fields.
    bool getIsInterfaceSpecificField() const {
      if (stagrLoc_ != "center") {
        return true;
      }
      const std::string ending = "_cold";
      if (longName_.size() > ending.size()) {
        // check if last `ending.size` characters of longName equal ending
        if (std::equal(longName_.begin() + longName_.size() - ending.size(),
                       longName_.end(), ending.begin())) {
          return true;
        }
      }
      return false;
    }

    // Set functions (strings)
    // -----------------------
    void setShrtName(std::string shrtName) {shrtName_ = shrtName;}
    void setDataKind(std::string dataKind) {dataKind_ = dataKind;}
    void setStagrLoc(std::string stagrLoc) {stagrLoc_ = stagrLoc;}
    void setMathSpac(std::string mathSpac) {mathSpac_ = mathSpac;}
    void setVarUnits(std::string varUnits) {varUnits_ = varUnits;}
    void setInOuName(std::string inOuName) {inOuName_ = inOuName;}
    void setInOuFile(std::string inOuFile) {inOuFile_ = inOuFile;}
    void setIntrpTyp(std::string intrpTyp) {intrpTyp_ = intrpTyp;}
    void setIntrpMsk(std::string intrpMsk) {intrpMsk_ = intrpMsk;}

    // Set number of levels
    // --------------------
    void setNumLevls(int numLevls) {numLevls_ = numLevls;}
    void setNumLevls(std::string numLevls) {
      if (numLevls == "full") {
        numLevls_ = nlev_;
      } else if (numLevls == "half") {
        numLevls_ = nlev_ + 1;
      } else if (numLevls == "halfplusone") {
        numLevls_ = nlev_ + 2;
      } else {
        try {
          numLevls_ = std::stoi(numLevls);
        } catch (std::invalid_argument& e) {
          ABORT("FieldMetadata::setFieldNumLevls levels neither full, half or an integer");
        }
      }
    }

    // Set tracer
    // ----------
    void setIsTracer(bool isTracer) {isTracer_ = isTracer;}
    void setIsTracer(std::string tracer) {
      if (tracer == "true") {
         isTracer_ = true;
      } else if (tracer == "false") {
         isTracer_ = false;
      } else {
        ABORT("FieldMetadata::setIsTracer tracer must be true or false");
      }
    }

    // Validity macro
    void validateVariable(std::vector<std::string> validOptions, std::string choice) const {
      auto result = std::find(validOptions.begin(), validOptions.end(), choice);
      if (result == std::end(validOptions)) {
        ABORT("FieldMetadata::validate For long name " + longName_ + " invalid kind: " + choice);
      }
    }

    // Check validity of choices
    void validate() const {
      this->validateVariable(ValidDataKind_, dataKind_);
      this->validateVariable(ValidIntrpTyp_, intrpTyp_);
      this->validateVariable(ValidMathSpac_, mathSpac_);
      this->validateVariable(ValidStagrLoc_, stagrLoc_);
    }

   private:
    // Picked up from default file
    std::string longName_;
    std::string shrtName_;
    std::string dataKind_;
    std::string stagrLoc_;
    int numLevls_;
    std::string mathSpac_;
    bool isTracer_;

    // Picked up from both default and override file
    std::string varUnits_;

    // Only picked up from override file
    std::string inOuName_;
    std::string inOuFile_;
    std::string intrpTyp_;
    std::string intrpMsk_;

    // Number of levels for the model
    int nlev_;

    // Valid choices
    const std::vector<std::string> ValidDataKind_ = {"double", "integer"};
    const std::vector<std::string> ValidIntrpTyp_ = {"integer", "nearest", "default"};
    const std::vector<std::string> ValidMathSpac_ = {"vector", "magnitude", "direction"};
    const std::vector<std::string> ValidStagrLoc_ = {"center", "eastwest", "northsouth", "corner"};

    // Print method
    void print(std::ostream & os) const {
      os << std::endl << "   Long name: " << longName_;
      os << std::endl << "   Short name: " << shrtName_;
      os << std::endl << "   Units: " << varUnits_;
      os << std::endl << "   Kind: " << dataKind_;
      os << std::endl << "   Horizontal stagger location: " << stagrLoc_;
      os << std::endl << "   Levels: " << numLevls_;
      os << std::endl << "   Space: " << mathSpac_;
      os << std::endl << "   Tracer: " << isTracer_;
      os << std::endl << "   IO name: " << inOuName_;
      os << std::endl << "   IO file: " << inOuFile_;
      os << std::endl << "   Interpolation type: " << intrpTyp_;
      os << std::endl << "   Interpolation source-point mask: " << intrpMsk_;
    }
  };

  // -----------------------------------------------------------------------------------------------

  class FieldsMetadata : public util::Printable {
   public:
    typedef FieldsMetadataParameters Parameters_;
    FieldsMetadata(const Parameters_ &, int &);

    // Get FieldMetadata from any of the potential field names
    FieldMetadata getFieldMetadata(const std::string &) const;

    // Get levels from any of the potential field names
    size_t getLevels(const std::string &) const;

    // Get long name from any of the potential field names
    oops::Variables getLongNameFromAnyName(const oops::Variables &) const;
    std::string getLongNameFromAnyName(const std::string &) const;

    // Filter out any fields that are specific to the fv3-jedi interface, returns a
    // Variables containing fields for passing to JEDI
    oops::Variables removeInterfaceSpecificFields(const oops::Variables &) const;

   private:
    std::map<std::string, FieldMetadata> fieldsMetadata_;

    // Print method
    void print(std::ostream & os) const {
      os << std::endl << " List of field meta data available: \n";
      for (const auto& ke : fieldsMetadata_) {
        os << std::endl << "  Key = " << ke.first << ":" << ke.second << "\n";
      }
    }
  };

  // -----------------------------------------------------------------------------------------------

}  // namespace fv3jedi
