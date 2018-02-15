#!/bin/sh --login

# script to build FMS from https://github.com/NOAA-GFDL/FMS on Theia
# 1) git clone https://github.com/NOAA-GFDL/FMS.git
# 2) cd FMS
# 3) run this script

source $MODULESHOME/init/sh

export CPPFLAGS="-Duse_netCDF -Duse_libMPI -DINTERNAL_FILE_NML -fPIC"
export PROJECT_DIR=$(groups | cut -d" " -f1)

module purge
module use -a /home/fms/local/modulefiles
module load fre/bronx-11
module load intel/15.6.233
module load netcdf/4.3.0
module load hdf5/1.8.14
module load impi/5.1.2.150
module list
list_paths -o pathnames_fms .
make clean
mkmf -m Makefile -p libfms.a -t $FRE_COMMANDS_HOME/site/$FRE_SYSTEM_SITE/intel.mk -c "$CPPFLAGS" -Iinclude -Impp/include pathnames_fms
make OPENMP=Y DEBUG=Y libfms.a

exit 0
