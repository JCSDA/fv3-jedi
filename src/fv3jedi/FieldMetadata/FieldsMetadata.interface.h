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
  void fields_metadata_get_field_f(const FieldsMetadata* fieldsMetadata, const char fieldIONameC[],
                                   char fieldNameC[], char kindC[], int& levels,
                                   char longNameC[], char SpaceC[], char staggerLocC[],
                                   bool& tracer, char unitsC[], char interpTypeC[],
                                   char io_fileC[]);
}

}  // namespace fv3jedi
