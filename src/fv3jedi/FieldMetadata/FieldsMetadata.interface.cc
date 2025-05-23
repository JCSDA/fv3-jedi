/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cstring>
#include <iostream>
#include <string>

#include "eckit/exception/Exceptions.h"

#include "oops/util/abor1_cpp.h"

#include "fv3jedi/FieldMetadata/FieldsMetadata.h"
#include "fv3jedi/FieldMetadata/FieldsMetadata.interface.h"

namespace fv3jedi {

  // -----------------------------------------------------------------------------------------------

  void checkStringLength(const std::string strCheck) {
    unsigned fortranStrLen = 2048;
    if (strCheck.size() >= fortranStrLen) {
      ABORT("FieldMetadataInterface.check_string_length "+strCheck+" exceeds length from Fortran");
    }
  }

  // -----------------------------------------------------------------------------------------------

  void get_field_metadata_f(const FieldsMetadata* fieldsMetadata,
                            const char longNameC[], char varUnitsC[], char dataKindC[],
                            bool& tracer, int & levels, char mathSpacC[]) {
    // Get meta data for requested field
    const std::string longName(longNameC);
    FieldMetadata fieldMetadata = fieldsMetadata->getFieldMetadata(longName);

    // Bool, int outputs
    levels = fieldMetadata.getNumLevls();
    tracer = fieldMetadata.getIsTracer();

    // Prepare char outputs
    std::string varUnits = fieldMetadata.getVarUnits();
    std::string dataKind = fieldMetadata.getDataKind();
    std::string mathSpac = fieldMetadata.getMathSpac();

    // Check string lengths
    checkStringLength(varUnits);
    checkStringLength(dataKind);
    checkStringLength(mathSpac);

    // Fill char outputs
    std::copy(varUnits.begin(), varUnits.end(), varUnitsC);
    std::copy(dataKind.begin(), dataKind.end(), dataKindC);
    std::copy(mathSpac.begin(), mathSpac.end(), mathSpacC);
  }

  // -----------------------------------------------------------------------------------------------

}  // namespace fv3jedi
