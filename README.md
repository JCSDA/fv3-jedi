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

