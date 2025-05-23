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

  FieldsMetadata::FieldsMetadata(const int nlev) : longNames_() {
    // Set the default metadata
    // ------------------------
    setMetadata(fieldsMetadata_, nlev);

    // Create vector of the field long names
    // -------------------------------------
    for (const auto & field : fieldsMetadata_) {
      longNames_.push_back(field.first);
    }
  }

  // -----------------------------------------------------------------------------------------------

  FieldMetadata FieldsMetadata::getFieldMetadata(const std::string & longName) const {
    // Check that fieldsMetadata_ has longName in the keys and abort if not
    ASSERT_MSG(fieldsMetadata_.find(longName) != fieldsMetadata_.end(),
               "FieldMetadata error. Field \"" + longName + "\" not found in map. Ensure that " +
               "the field is listed in FieldMetadataDefault.h");
    // Return Field Metadata
    return fieldsMetadata_.find(longName)->second;
  }

  // -----------------------------------------------------------------------------------------------

  size_t FieldsMetadata::getLevels(const std::string & longName) const {
    // Get the element
    const FieldMetadata field = this->getFieldMetadata(longName);
    // Return number of levels
    return field.getNumLevls();
  }

  // -----------------------------------------------------------------------------------------------

}  // namespace fv3jedi

// -------------------------------------------------------------------------------------------------
