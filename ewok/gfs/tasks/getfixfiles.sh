#!/bin/bash

# (C) Copyright 2020-2021 UCAR
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
cp ${FV3REPO}/test/Data/fv3files/field_table_gfdl ${FIXDIR}/field_table_gfdl
cp ${FV3REPO}/test/Data/fieldmetadata/gfs-restart.yaml ${FIXDIR}/gfs-restart.yaml

cp ${FV3REPO}/test/Data/fv3files/input_gfs_${RESOL}_${LAYOUT}.nml ${FIXDIR}/model_gfs_${RESOL}.nml

