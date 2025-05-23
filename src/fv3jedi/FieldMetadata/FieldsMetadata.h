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

#include "oops/util/abor1_cpp.h"
#include "oops/util/Printable.h"

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
    std::string getDataKind() const {return dataKind_;}
    std::string getMathSpac() const {return mathSpac_;}
    std::string getVarUnits() const {return varUnits_;}

    // Set functions (strings)
    // -----------------------
    void setDataKind(std::string dataKind) {dataKind_ = dataKind;}
    void setMathSpac(std::string mathSpac) {mathSpac_ = mathSpac;}
    void setVarUnits(std::string varUnits) {varUnits_ = varUnits;}

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
      this->validateVariable(ValidMathSpac_, mathSpac_);
    }

   private:
    // Picked up from default file
    std::string longName_;
    std::string dataKind_;
    int numLevls_;
    std::string mathSpac_;
    bool isTracer_;

    // Picked up from both default and override file
    std::string varUnits_;

    // Number of levels for the model
    int nlev_;

    // Valid choices
    const std::vector<std::string> ValidDataKind_ = {"double", "integer"};
    const std::vector<std::string> ValidMathSpac_ = {"vector", "magnitude", "direction"};

    // Print method
    void print(std::ostream & os) const {
      os << std::endl << "   Long name: " << longName_;
      os << std::endl << "   Units: " << varUnits_;
      os << std::endl << "   Kind: " << dataKind_;
      os << std::endl << "   Levels: " << numLevls_;
      os << std::endl << "   Space: " << mathSpac_;
      os << std::endl << "   Tracer: " << isTracer_;
    }
  };

  // -----------------------------------------------------------------------------------------------

  class FieldsMetadata : public util::Printable {
   public:
    explicit FieldsMetadata(const int);

    // Get FieldMetadata from any of the potential field names
    FieldMetadata getFieldMetadata(const std::string &) const;

    // Get levels from any of the potential field names
    size_t getLevels(const std::string &) const;

    // Function to return all the long names
    const std::vector<std::string> & getLongNames() const {return longNames_;}

   private:
    std::map<std::string, FieldMetadata> fieldsMetadata_;
    std::vector<std::string> longNames_;

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
