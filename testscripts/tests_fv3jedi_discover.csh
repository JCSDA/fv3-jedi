#!/bin/csh -fx
# ------------------------------
#SBATCH -A g0613
#SBATCH --export=NONE
#SBATCH --job-name=jedifv3tests
#SBATCH --output=jedifv3tests.log.o%j
#SBATCH --nodes=1
#SBATCH --qos=debug
#SBATCH --tasks-per-node=6
#SBATCH --constraint=hasw
#SBATCH --time=01:00:00

#Script for testing fv3/jedi on the NASA Center for Climate Similation (NCCS) Discover cluster

#Shell
source /usr/share/modules/init/csh

#Jedi install config
setenv JEDI_ROOT /discover/nobackup/drholdaw/Jedi/

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

   #Boost/Eigen include dirs
   setenv BOOST_ROOT /discover/nobackup/drholdaw/Jedi/Jedi_Shared/boost/1.66.0_int/include/
   setenv EIGEN3_PATH /discover/nobackup/drholdaw/Jedi/Jedi_Shared/eigen/3.3.4_int/include/

else if ($1 == "GCC" || $1 == "gcc" || $1 == "GNU" || $1 == "gnu") then
   
   #Modules GCC
   module purge
   module load other/cmake-3.8.2
   module load other/comp/gcc-7.2
   module load other/mpi/openmpi/3.0.0-gcc-7.2
   setenv BASEDIR /discover/swdev/mathomp4/Baselibs/ESMA-Baselibs-4.0.10/x86_64-unknown-linux-gnu/gfortran_7.2.0-openmpi_3.0.0/Linux
   setenv MPIEXEC `which mpirun`
   setenv CPCcomp mpicxx
   setenv CCcomp mpicc
   setenv F90comp mpifort

   module load boost/1.66.0_gcc
   module load eigen/3.3.4_gcc

   setenv JEDI_BUILD $JEDI_ROOT/build_gcc

   #Boost/Eigen include dirs
   setenv BOOST_ROOT /discover/nobackup/drholdaw/Jedi/Jedi_Shared/boost/1.66.0_gcc/include/
   setenv EIGEN3_PATH /discover/nobackup/drholdaw/Jedi/Jedi_Shared/eigen/3.3.4_gcc/include/

else
   
   echo "No compiler specified"
   exit()

endif

#Run the tests
cd $JEDI_BUILD/fv3-jedi/test/

ctest -VV
