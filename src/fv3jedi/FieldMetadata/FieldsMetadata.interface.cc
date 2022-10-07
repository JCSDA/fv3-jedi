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

  void fields_metadata_get_field_f(const FieldsMetadata* fieldsMetadata,
                                   const char longshortioNameC[], char longNameC[],
                                   char shrtNameC[], char varUnitsC[], char dataKindC[],
                                   bool& tracer, char stagrLocC[], int & levels, char mathSpacC[],
                                   char inOuNameC[], char inOuFileC[], char intrpTypC[],
                                   char intrpMskC[]) {
    // Get meta data for requested field
    const std::string longshortioName(longshortioNameC);
    FieldMetadata fieldMetadata = fieldsMetadata->getField(longshortioName);

    // Bool, int outputs
    levels = fieldMetadata.getNumLevls();
    tracer = fieldMetadata.getIsTracer();

    // Prepare char outputs
    std::string longName = fieldMetadata.getLongName();
    std::string shrtName = fieldMetadata.getShrtName();
    std::string varUnits = fieldMetadata.getVarUnits();
    std::string dataKind = fieldMetadata.getDataKind();
    std::string stagrLoc = fieldMetadata.getStagrLoc();
    std::string mathSpac = fieldMetadata.getMathSpac();
    std::string inOuName = fieldMetadata.getInOuName();
    std::string inOuFile = fieldMetadata.getInOuFile();
    std::string intrpTyp = fieldMetadata.getIntrpTyp();
    std::string intrpMsk = fieldMetadata.getIntrpMsk();

    // Check string lengths
    checkStringLength(longName);
    checkStringLength(shrtName);
    checkStringLength(varUnits);
    checkStringLength(dataKind);
    checkStringLength(stagrLoc);
    checkStringLength(mathSpac);
    checkStringLength(inOuName);
    checkStringLength(inOuFile);
    checkStringLength(intrpTyp);
    checkStringLength(intrpMsk);

    // Fill char outputs
    std::copy(longName.begin(), longName.end(), longNameC);
    std::copy(shrtName.begin(), shrtName.end(), shrtNameC);
    std::copy(varUnits.begin(), varUnits.end(), varUnitsC);
    std::copy(dataKind.begin(), dataKind.end(), dataKindC);
    std::copy(stagrLoc.begin(), stagrLoc.end(), stagrLocC);
    std::copy(mathSpac.begin(), mathSpac.end(), mathSpacC);
    std::copy(inOuName.begin(), inOuName.end(), inOuNameC);
    std::copy(inOuFile.begin(), inOuFile.end(), inOuFileC);
    std::copy(intrpTyp.begin(), intrpTyp.end(), intrpTypC);
    std::copy(intrpMsk.begin(), intrpMsk.end(), intrpMskC);
  }

}  // namespace fv3jedi
