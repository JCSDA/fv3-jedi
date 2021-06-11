#!/bin/bash

# (C) Copyright 2017-2020 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

# remove the output file, if it exists
rm -f testoutput/${COMPARE_TESTNAME}.run

# run the executable
echo ""
echo "==============================================================================="
echo "Running test executable"
echo "Command: "
echo ${MPI_CMD} $1 $2 testoutput/${COMPARE_TESTNAME}.run
echo " "
echo "==============================================================================="
${MPI_CMD} $1 $2 testoutput/${COMPARE_TESTNAME}.run
e=$?
if [[ $e -gt 0 ]]; then
    echo -e "Failed to run executable. Error code: $e \n"
    exit $e
fi

# run compare, if needed
if [[ $SKIP_COMPARE == "FALSE" ]]; then
    echo ""
    echo "==============================================================================="
    echo "Running compare script"
    echo "Command: "
    echo $COMPARE_SCRIPT testoutput/${COMPARE_TESTNAME}.run testoutput/${COMPARE_REFERENCE}.ref ${COMPARE_TOL_F} ${COMPARE_TOL_I}
    echo " "
    echo "==============================================================================="

    $COMPARE_SCRIPT testoutput/${COMPARE_TESTNAME}.run testoutput/${COMPARE_REFERENCE}.ref ${COMPARE_TOL_F} ${COMPARE_TOL_I}
    e=$?
    if [[ $e -gt 0 ]]; then
        echo -e "Failed in the COMPARE step. Error code: $e \n"
        exit $e
    fi
    echo -e "PASSED \n"
fi
