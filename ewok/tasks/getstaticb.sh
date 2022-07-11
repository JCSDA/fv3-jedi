#!/bin/bash

# (C) Copyright 2020-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set -eux

echo "FV3DATA = ${FV3DATA}"
echo "STATBDIR = ${STATBDIR}"
echo "BALFILE = ${BALFILE}"
echo "PREFIX2D = ${PREFIX2D}"
echo "PREFIX3D = ${PREFIX3D}"

mkdir -p ${STATBDIR}

cp ${FV3DATA}/inputs/nmcbalance/${BALFILE}  ${STATBDIR}

cp ${FV3DATA}/bump/${PREFIX2D}*.nc  ${STATBDIR}
cp ${FV3DATA}/bump/${PREFIX3D}*.nc  ${STATBDIR}

