/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <string>

#include "oops/util/abor1_cpp.h"

#include "fv3jedi/FieldMetadata/FieldsMetadata.h"
#include "fv3jedi/FieldMetadata/FieldsMetadata.interface.h"

namespace fv3jedi {

  void checkStringLength(const std::string strCheck) {
    unsigned fortranStrLen = 2048;
    if (strCheck.size() >= fortranStrLen) {
      ABORT("FieldMetadataInterface.check_string_length "+strCheck+" exceeds length from Fortran");
    }
  }

  void fields_metadata_get_field_f(const FieldsMetadata* fieldsMetadata, const char fieldIONameC[],
                                   char fieldNameC[], char kindC[], int& levels,
                                   char longNameC[], char spaceC[], char staggerLocC[],
                                   bool& tracer, char unitsC[], char interpTypeC[],
                                   char io_fileC[]) {
    // Get meta data for requested field
    oops::Log::trace() << "Calling FieldsMetaData.GetField for " << fieldIONameC << std::endl;
    const std::string fieldIOName(fieldIONameC);
    FieldMetadata fieldMetadata = fieldsMetadata->getField(fieldIOName);

    // Bool, int outputs
    levels = fieldMetadata.getLevels();
    tracer = fieldMetadata.getTracer();

    // Prepare char outputs
    std::string fieldName = fieldMetadata.getFieldName();
    std::string kind = fieldMetadata.getKind();
    std::string longName = fieldMetadata.getLongName();
    std::string space = fieldMetadata.getSpace();
    std::string staggerLoc = fieldMetadata.getStaggerLoc();
    std::string units = fieldMetadata.getUnits();
    std::string interpType = fieldMetadata.getInterpType();
    std::string io_file = fieldMetadata.getIOFile();

    // Check string lengths
    checkStringLength(fieldName);
    checkStringLength(kind);
    checkStringLength(longName);
    checkStringLength(space);
    checkStringLength(staggerLoc);
    checkStringLength(units);
    checkStringLength(interpType);
    checkStringLength(io_file);

    // Fill char outputs
    std::copy(fieldName.begin(), fieldName.end(), fieldNameC);
    std::copy(kind.begin(), kind.end(), kindC);
    std::copy(longName.begin(), longName.end(), longNameC);
    std::copy(space.begin(), space.end(), spaceC);
    std::copy(staggerLoc.begin(), staggerLoc.end(), staggerLocC);
    std::copy(units.begin(), units.end(), unitsC);
    std::copy(interpType.begin(), interpType.end(), interpTypeC);
    std::copy(io_file.begin(), io_file.end(), io_fileC);
  }

}  // namespace fv3jedi
