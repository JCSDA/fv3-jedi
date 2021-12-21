#!/bin/bash

# (C) Copyright 2020-2021 UCAR
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set -eux

echo "FV3REPO = ${FV3REPO}"
echo "FIXDIR  = ${FIXDIR}"

cp ${FV3REPO}/test/Data/fv3files/akbk64.nc4       ${FIXDIR}/akbk64.nc4
cp ${FV3REPO}/test/Data/fv3files/fmsmpp.nml       ${FIXDIR}/fmsmpp.nml
cp ${FV3REPO}/test/Data/fv3files/field_table_gfdl ${FIXDIR}/field_table_gfdl
cp ${FV3REPO}/test/Data/fieldsets/dynamics.yaml   ${FIXDIR}/dynamics.yaml
cp ${FV3REPO}/test/Data/fieldsets/ufo.yaml        ${FIXDIR}/ufo.yaml

echo "RESOL = ${RESOL}"

cp ${FV3REPO}/test/Data/fv3files/input_gfs_${RESOL}.nml ${FIXDIR}/model_gfs_${RESOL}.nml

