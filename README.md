[![AWS-gnu](https://codebuild.us-east-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoiUEFLOXlGNjI2OUYybE9RTncyWS9hRjhXMGVXOXdKTUh1aVNuempFRHRMRjA1dmd5UDRienNaa2x6MDFHcWhhYmwvUmtRcnA5Y0FsbG8zS0xib3NtYTkwPSIsIml2UGFyYW1ldGVyU3BlYyI6Iks4Zkx2S3FKZGhwbXZzeEwiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=develop)](https://us-east-1.console.aws.amazon.com/codesuite/codebuild/projects/automated-testing-fv3-gnu/history)
[![AWS-intel](https://codebuild.us-east-1.amazonaws.com/badges?uuid=eyJlbmNyeXB0ZWREYXRhIjoiVHRMaWZsV1VSd2luTzM1eEhwS0VjRWx3ajNHSTBkRThLZzFaWnJVQ3VBTkJGY0wwazllSHJXVVRVM3BLTlo5YUtWZ0N5Z2hNWTlOU1M1WWJMVklBeEZNPSIsIml2UGFyYW1ldGVyU3BlYyI6IlFpOVlIWFFGOXo3UlhnQlAiLCJtYXRlcmlhbFNldFNlcmlhbCI6MX0%3D&branch=develop)](https://us-east-1.console.aws.amazon.com/codesuite/codebuild/projects/automated-testing-fv3-intel/history)

FV3JEDI interfaces for OOPS

(C) Copyright 2017 UCAR.

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

--- Requirements ---

See OOPS requirements

--- Building ---

The variables ${SRC_OOPS}, ${SRC_MODEL} and ${BUILD} below must be defined for your
environement.
Note: It is good practice to build the code outside of the source tree.

The lines below can be copied into a script or executed manually:

Define environment

    export SRC_OOPS=/path/to/source/oops
    export SRC_MODEL=/path/to/source/fv3jedi
    export BUILD=/path/to/build

    export PATH=${PATH}:${SRC_OOPS}/ecbuild/bin

Build OOPS first

    rm -rf ${BUILD}/oops; mkdir ${BUILD}/oops; cd ${BUILD}/oops
    ecbuild --build=release ${SRC_OOPS}
    make -j4

Model specific build

    rm -rf ${BUILD}/fv3jedi; mkdir ${BUILD}/fv3jedi; cd ${BUILD}/fv3jedi
    ecbuild -DOOPS_PATH=${BUILD}/oops --build=release ${SRC_MODEL}
    make -j4

For testing the build:

    cd ${BUILD}/fv3jedi
    ctest

--- Working with OOPS ---

After the code has been built successfully once, it is enough to re-run the make
command only for re-compiling the code after modifications.

