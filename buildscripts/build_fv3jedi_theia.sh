#!/bin/sh --login

# run as qsub -d $PWD -l nodes=1:ppn=24 -l walltime=00:30:00 -q debug -a fv3-cpu -j oe -o oops_fv3jedi.out build_fv3jedi_theia.sh

# Build OOPS or skip it.
SKIP_BUILD_OOPS="NO"

SRC_ROOT="$(pwd)/.."
SRC_OOPS=$SRC_ROOT/jedi
SRC_MODEL=$SRC_ROOT/fv3-jedi
SRC_FMS=$SRC_ROOT/FMS

BUILD=$SRC_ROOT/build

# Build OOPS (or skip it)
if [ ${SKIP_BUILD_OOPS:-"NO"} = "NO" ]; then

    sh $SRC_OOPS/build_oops_theia.sh
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "OOPS failed building or some tests failed"
        exit $rc
    fi

else

    echo "Skipping OOPS build"
    echo "Expecting OOPS built under:"
    echo "$BUILD/oops"
    if [ "$(ls -A $BUILD/oops)" ]; then
        echo "ls $BUILD/oops"
        ls $BUILD/oops
    else
        echo "$BUILD/oops is empty or does not exist, ABORT!"
        exit -1
    fi

fi

# Clean modules
source $MODULESHOME/init/sh

module purge

# Load Intel compilers, NetCDF and HDF5 libraries
module load intel/15.6.233
module load impi/5.1.2.150
module load hdf5/1.8.14 netcdf/4.3.0

# Load cmake
module use -a /contrib/modulefiles
module load cmake

# Load Boost and Eigen
module use -a /contrib/da/modulefiles
module load boost
module load eigen

module list

# Need eigen3 and netcdf libraries
export EIGEN3_INCLUDE_DIR=$EIGEN_ROOT
export NETCDF_LIBRARIES="${NETCDF}/lib/libnetcdf.a;${NETCDF}/lib/libnetcdff.a"

# Need correct MPIEXEC on Theia
MPIEXEC=$(which mpirun)

# Add ecbuild/bin to path
export PATH=${PATH}:${SRC_OOPS}/ecbuild/bin

export FMS_LIBRARIES="$SRC_FMS/libfms.a"
export FMS_INCLUDES="$SRC_FMS;$SRC_FMS/include"

# Build FV3JEDI

rm -rf $BUILD/fv3-jedi; mkdir -p $BUILD/fv3-jedi; cd $BUILD/fv3-jedi

ecbuild \
    --build=release \
    -DCMAKE_CXX_COMPILER=mpiicpc \
    -DCMAKE_C_COMPILER=mpiicc \
    -DCMAKE_Fortran_COMPILER=mpiifort \
    -DBOOST_ROOT=$BOOST_ROOT -DBoost_NO_SYSTEM_PATHS=ON \
    -DOOPS_PATH=$BUILD/oops \
    -DMPIEXEC=$MPIEXEC \
    $SRC_MODEL

make -j2

# For testing the build:
/bin/cp test/testinput/input.nml test # input.nml needs to be one level up
/bin/cp -R $SRC_MODEL/test/fv3jedi_geom/INPUT test
mkdir test/RESTART
ctest -VV

exit 0
