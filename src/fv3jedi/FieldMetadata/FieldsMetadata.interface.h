/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "fv3jedi/FieldMetadata/FieldsMetadata.h"

namespace fv3jedi {

extern "C" {
  void fields_metadata_get_field_f(const FieldsMetadata* fieldsMetadata,
                                   const char longshortioNameC[], char longNameC[],
                                   char shrtNameC[], char varUnitsC[], char dataKindC[],
                                   bool& tracer, char stagrLocC[], int & levels, char mathSpacC[],
                                   char inOuNameC[], char inOuFileC[], char intrpTypC[]);

}

}  // namespace fv3jedi
