/*
 * (C) Copyright 2020-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

// -------------------------------------------------------------------------------------------------

#include <iostream>
#include <map>
#include <utility>

#include "eckit/exception/Exceptions.h"
#include "fv3jedi/FieldMetadata/FieldsMetadata.h"
#include "fv3jedi/FieldMetadata/FieldsMetadataDefault.h"

// -------------------------------------------------------------------------------------------------

namespace fv3jedi {

  // -----------------------------------------------------------------------------------------------

  FieldsMetadata::FieldsMetadata(const Parameters_ & params, int & nlev) {
    // Set the default metadata
    // ------------------------
    setDefaults(fieldsMetadata_, nlev);

    // If necessary open the override file and replace optionally replaceable metadata
    // -------------------------------------------------------------------------------
    if (params.override.value() != boost::none) {
      // Open override Yaml file
      eckit::PathName overrideMetadataFile(*params.override.value());
      const eckit::YAMLConfiguration overrideMetadataDict(overrideMetadataFile);

      // Create parameters object from the YAML file
      FieldsOverrideParameters fieldsOverrideParamsDict;
      fieldsOverrideParamsDict.validateAndDeserialize(overrideMetadataDict);

      // List of fields
      const std::vector<FieldOverrideParameters> fieldsOverride =
                                                            fieldsOverrideParamsDict.fields.value();

      // Loop over list within field metadata key (main list of override field metadata)
      for (const auto & fieldOverride : fieldsOverride) {
        // Get the long name from the parameters
        const std::string longName = fieldOverride.longName.value();

        // Check that key is in the map
        ASSERT_MSG(fieldsMetadata_.find(longName) != fieldsMetadata_.end(),
                   "FieldMetadata::FieldsMetadata: Trying to override " + longName +
                   " but this long name does not exist in the metadata.");

        // Get pointer to object
        FieldMetadata& fieldMetadata = fieldsMetadata_.find(longName)->second;

        // Units
        if (fieldOverride.varUnits.value() != boost::none) {
          fieldMetadata.setVarUnits(*fieldOverride.varUnits.value());
        }

        // IO name
        if (fieldOverride.InOuName.value() != boost::none) {
          fieldMetadata.setInOuName(*fieldOverride.InOuName.value());
        }

        // IO file
        if (fieldOverride.InOuFile.value() != boost::none) {
          fieldMetadata.setInOuFile(*fieldOverride.InOuFile.value());
        }

        // Interpolation type
        if (fieldOverride.IntrpTyp.value() != boost::none) {
          fieldMetadata.setIntrpTyp(*fieldOverride.IntrpTyp.value());
        }

        // Interpolation source-point mask
        if (fieldOverride.IntrpMsk.value() != boost::none) {
          fieldMetadata.setIntrpMsk(*fieldOverride.IntrpMsk.value());
        }
      }
    }

    // Check for duplicated short/io name
    // ----------------------------------
    // Get a vector will all short/io names
    std::vector<std::string> allShrtNames;
    std::vector<std::string> allInOuNames;
    for (const auto& ke : fieldsMetadata_) {
      allShrtNames.push_back(ke.second.getShrtName());
      allInOuNames.push_back(ke.second.getInOuName());
    }

    // Sort the vectors
    sort(allShrtNames.begin(), allShrtNames.end());
    sort(allInOuNames.begin(), allInOuNames.end());

    // Identify if there is a duplicate in either
    const auto duplicateShrt = std::adjacent_find(allShrtNames.begin(), allShrtNames.end());
    const auto duplicateInOu = std::adjacent_find(allInOuNames.begin(), allInOuNames.end());

    // Report duplicates and fail if needed
    if (duplicateShrt != allShrtNames.end()) {
      ABORT("FieldMetadata::FieldsMetadata: Short name "+*duplicateShrt+" is duplicated.");
    }
    if (duplicateInOu != allInOuNames.end()) {
      ABORT("FieldMetadata::FieldsMetadata: IO name "+*duplicateInOu+" is duplicated.");
    }
  }

  // -----------------------------------------------------------------------------------------------

  FieldMetadata FieldsMetadata::getField(const std::string & longshortio) const {
    // Get longname in case incoming is not already long name
    const std::string longName = this->getLongNameFromAnyName(longshortio);
    // Return Field Metadata
    return fieldsMetadata_.find(longName)->second;
  }

  // -----------------------------------------------------------------------------------------------

  size_t FieldsMetadata::getLevels(const std::string & longshortio) const {
    // Get the element
    const FieldMetadata field = this->getField(longshortio);
    // Return number of levels
    return field.getNumLevls();
  }

  // -----------------------------------------------------------------------------------------------

  std::string FieldsMetadata::getLongNameFromAnyName(const std::string & longshortio) const {
    // Loop through the map and look for longshortio name
    for (const auto& ke : fieldsMetadata_) {
      if (longshortio == ke.second.getInOuName() || longshortio == ke.second.getShrtName() ||
          longshortio == ke.first) {
        return ke.second.getLongName();
      }
    }

    // Fail if return not hit
    ABORT("FieldMetadata::getLongNameFromAnyName: Searching for a field called "+longshortio+
          " in the long, short and io names but not found anywhere.");
    return "Failed";  // Avoids compilation warning
  }

  // -----------------------------------------------------------------------------------------------

  oops::Variables FieldsMetadata::getLongNameFromAnyName(const oops::Variables & vars) const {
    // Method takes as input an oops::Variables object containing either long names, short names or
    // io names. It returns oop::Variables objects with long names.

    // Convert variables to vector of strings
    const std::vector<std::string>& varsVec = vars.variables();

    // Vector of long name strings
    std::vector<std::string> longNameVec;

    // Iterate over vars and find equivalent long name
    for (auto &var : varsVec) {
      longNameVec.push_back(this->getLongNameFromAnyName(var));
    }

    // Return long name vars
    return oops::Variables(longNameVec);
  }

  // -----------------------------------------------------------------------------------------------

}  // namespace fv3jedi

// -------------------------------------------------------------------------------------------------
