#!/bin/bash

# (C) Copyright 2020-2022 UCAR
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set -eux

echo "FV3REPO = ${FV3REPO}"
echo "FIXDIR  = ${FIXDIR}"

echo "RESOL = ${RESOL}"
echo "NLEVS = ${NLEVS}"
echo "LAYOUT = ${LAYOUT}"

cp ${FV3REPO}/test/Data/fv3files/akbk${NLEVS}.nc4 ${FIXDIR}/akbk${NLEVS}.nc4
cp ${FV3REPO}/test/Data/fv3files/fmsmpp.nml       ${FIXDIR}/fmsmpp.nml
cp ${FV3REPO}/test/Data/fv3files/field_table_gmao ${FIXDIR}/field_table_gmao
cp ${FV3REPO}/test/Data/fieldmetadata/geos_cf.yaml ${FIXDIR}/geos_cf.yaml

# we are not running the model in fc or inner loop yet...
cp ${FV3REPO}/test/Data/fv3files/input_geos_${RESOL}_z72.nml ${FIXDIR}/model_geos_${RESOL}.nml

