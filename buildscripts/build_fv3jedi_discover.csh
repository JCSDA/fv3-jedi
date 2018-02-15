#!/bin/csh -f

#Script for compiling fv3/jedi on the NASA Center for Climate Similation (NCCS) Discover cluster

#Clear path
set path = /bin:/usr/bin

#Shell
source /usr/share/modules/init/csh

#Jedi install config
setenv JEDI_ROOT /discover/nobackup/drholdaw/Jedi/
setenv JEDI_SRC $JEDI_ROOT/jedi-bundle/

module purge
module load other/cmake-3.8.2
module use -a /discover/nobackup/drholdaw/Jedi/Jedi_Shared/modulefiles

if ($1 == "INT" || $1 == "Int" || $1 == "Intel"  || $1 == "intel") then

   module load other/comp/gcc-6.2
   module load comp/intel-18.0.1.163
   module load mpi/impi-18.0.1.163
   module load lib/mkl-18.0.1.163
   setenv BASEDIR /discover/swdev/mathomp4/Baselibs/ESMA-Baselibs-5.0.9/x86_64-unknown-linux-gnu/ifort_18.0.1.163-intelmpi_18.0.1.163/Linux/ 

   setenv MPIEXEC `which mpirun`
   setenv CPCcomp mpiicpc
   setenv CCcomp mpiicc
   setenv F90comp mpiifort
   
   module load boost/1.66.0_int
   module load eigen/3.3.4_int

   setenv JEDI_BUILD $JEDI_ROOT/build_int
   setenv GFDL_BUILD $JEDI_ROOT/build_gfdl_int

   #Boost/Eigen include dirs
   setenv BOOST_ROOT /discover/nobackup/drholdaw/Jedi/Jedi_Shared/boost/1.66.0_int/include/
   setenv EIGEN3_PATH /discover/nobackup/drholdaw/Jedi/Jedi_Shared/eigen/3.3.4_int/include/

else if ($1 == "GCC" || $1 == "gcc" || $1 == "GNU" || $1 == "gnu") then
   
   #Modules GCC
   module load other/comp/gcc-7.2  
   module load other/mpi/openmpi/3.0.0-gcc-7.2
   module load lib/mkl-18.0.1.163
   setenv BASEDIR /discover/swdev/mathomp4/Baselibs/ESMA-Baselibs-4.0.10/x86_64-unknown-linux-gnu/gfortran_7.2.0-openmpi_3.0.0/Linux
   setenv MPIEXEC `which mpirun`
   setenv CPCcomp mpicxx
   setenv CCcomp mpicc
   setenv F90comp mpifort

   module load boost/1.66.0_gcc
   module load eigen/3.3.4_gcc

   setenv JEDI_BUILD $JEDI_ROOT/build_gcc
   setenv GFDL_BUILD $JEDI_ROOT/build_gfdl_gcc

   #Boost/Eigen include dirs
   setenv BOOST_ROOT /discover/nobackup/drholdaw/Jedi/Jedi_Shared/boost/1.66.0_gcc/include/
   setenv EIGEN3_PATH /discover/nobackup/drholdaw/Jedi/Jedi_Shared/eigen/3.3.4_gcc/include/

else
   
   echo "No compiler specified"
   exit()

endif

#Source to be built
setenv SRC $JEDI_ROOT/fv3-jedi/

#NetCDF lib search path
setenv NETCDF $JEDI_BUILD/netcdf
setenv NETCDF_INCLUDE_DIRS $JEDI_BUILD/netcdf/include #Need in order not to auto redefine NETCDF_LIBRARIES
setenv NETCDF_LIBRARIES "$JEDI_BUILD/netcdf/lib/libnetcdf.a;$JEDI_BUILD/netcdf/lib/libnetcdff.a;${BASEDIR}/lib/libhdf5_hl.a;${BASEDIR}/lib/libhdf5.a;${BASEDIR}/lib/libcurl.a;/usr/lib64/libcrypto.so;/usr/lib64/libssl.so;${BASEDIR}/lib/libmfhdf.a;${BASEDIR}/lib/libdf.a;${BASEDIR}/lib/libjpeg.a"

#Set defs from GMAO builds
setenv COMPDEFS "-DTLADPRES"

#Add build to path
set path = (${path} ${JEDI_ROOT}/ecbuild/bin)

cd $JEDI_BUILD

if ($2 == "clean" || ! -d fv3-jedi/) then

   #Prepare build location
   rm -rf fv3-jedi
   mkdir -p fv3-jedi
   cd fv3-jedi
   
   #Prepare to make
   ecbuild \
       --build=debug \
       -DCMAKE_CXX_COMPILER=$CPCcomp \
       -DCMAKE_C_COMPILER=$CCcomp \
       -DCMAKE_Fortran_COMPILER=$F90comp \
       -DBOOST_ROOT=$BOOST_ROOT -DBoost_NO_SYSTEM_PATHS=ON \
       -DOOPS_PATH=$JEDI_BUILD/oops \
       -DNETCDF_INCLUDE_DIRS=$NETCDF_INCLUDE_DIRS \
       -DNETCDF_LIBRARIES=$NETCDF_LIBRARIES \
       -DNETCDF_PATH=$NETCDF \
       -DFMS_PATH=$GFDL_BUILD/fms/ \
       -DFV3_PATH=$GFDL_BUILD/fv3/ \
       -DUFO_PATH=$JEDI_BUILD/ufo \
       -DMPIEXEC=$MPIEXEC \
       -DCOMPDEFS=$COMPDEFS \
       $SRC

endif

#Build the model
cd $JEDI_BUILD/fv3-jedi
make VERBOSE=YES -j2

#Run the tests
#cd $JEDI_BUILD/fv3-jedi
#ctest
