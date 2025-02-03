#!/usr/bin/env bash
# this script sets up inputs and sets dummy variables so that the first part of the ufs-weather-model run_test.py script will
# work on a generic system and generate the experiment directory that is used to run the ufs forecast regression tests.
set -x

mkdir -p Data/ModelRunDirs/c48_001
mkdir -p Data/ModelRunDirs/c48_002
cp  Data/ModelRunDirs/UFS_warmstart_1/input.nml  Data/ModelRunDirs/c48_001

echo "&ensemble_nml
 ensemble_size = 2
/
" >> Data/ModelRunDirs/c48_001/input.nml
cp  Data/ModelRunDirs/c48_001/input.nml  Data/ModelRunDirs/c48_002
cd Data/ModelRunDirs/c48_001
mkdir INPUT
ln -s ../UFS_warmstart_1/* .
cd INPUT
cp -r ../../UFS_warmstart_1/INPUT/* .
for file in *tile*.nc; do new=`echo $file | sed 's/tile/ens_01.tile/g'`; ln -s $file $new; done
ln -s fv_core.res.nc fv_core.res.ens_01.nc

cd ../../c48_002
mkdir INPUT
ln -s ../UFS_warmstart_2/* .
cd INPUT
cp -r ../../UFS_warmstart_2/INPUT/* .
for file in *tile*.nc; do new=`echo $file | sed 's/tile/ens_02.tile/g'`; ln -s $file $new; done
ln -s fv_core.res.nc fv_core.res.ens_02.nc

