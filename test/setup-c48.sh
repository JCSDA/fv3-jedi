#!/usr/bin/env bash
# this script sets up inputs and sets dummy variables so that the first part of the ufs-weather-model run_test.py script will
# work on a generic system and generate the experiment directory that is used to run the ufs forecast regression tests.
set -x

# create the input data directory which will be used by run_test.py
mkdir input-data
export INPUTDATA_ROOT=$PWD/input-data
export INPUTDATA_ROOT_WW3=$PWD/input-data
export INPUTDATA_ROOT_BMIC=$PWD/input-data

# script command line variables pointing to the build directory of ufs-bundle, the source directory of ufs-weather-model
# and the location of the fv3jedi data director
export builddir=$1
export ufsdir=$2
export fv3jedidata=$3

# Get date used by latest regression tests of the wm for aws sync below
baseline=`grep baseline $ufsdir/tests/logs/RegressionTests_hera.log | grep "\/control_c48_intel" | awk -F "/" '{print $7}'`

cd input-data
# pull the latest fixe files, input data and restart files from aws backup of rt.sh
aws s3 sync --no-sign-request s3://noaa-ufs-regtests-pds/input-data-20221101/FV3_fix FV3_fix
aws s3 sync --no-sign-request s3://noaa-ufs-regtests-pds/input-data-20221101/FV3_input_data48 FV3_input_data48
aws s3 sync --no-sign-request s3://noaa-ufs-regtests-pds/$baseline/control_c48_intel/RESTART RESTART

# set dummy values here so that run_test.py will think it is running on a supported platform
export RT_SUFFIX=
export BL_SUFFIX=
export SCHEDULER=slurm
export ACCNR=da-cpu
export QUEUE=batch
export PARTITION=
export ROCOTO=false
export ECFLOW=false
export REGRESSIONTEST_LOG=$ufsdir/rlog
export LOG_DIR=$ufsdir/log
export DEP_RUN=
export skip_check_results=false
export delete_rundir=false
export WLCLK=30
export JOB_NR="001"
touch fv3_001.exe
touch modules.fv3_001.lua
export PATHTR=$ufsdir
export RT_COMPILER=intel
export MACHINE_ID=hera
mkdir -p $LOG_DIR
touch $LOGDIR/job_001_timestamp.txt
cd $ufsdir/tests
touch fv3_001.exe
touch modules.fv3_001.lua
source vars

if [ "$(uname)" == "Darwin" ]; then
  hash gsed 2>/dev/null || { echo >&2 "GNU SED (gsed) required on macOS, but not installed. Aborting."; exit 1; }
  SED=gsed
else
  SED=sed
fi

# short circuit the job submission here. We only care about getting the input data set up
${SED} -i "/submit_and_wait/a exit 0" rt_utils.sh

# Generate the control_c48_intel regression test experiment directory for a cold start
./run_test.sh $PWD $fv3jedidata control_c48 001 001

cd $fv3jedidata/control_c48_intel

# Now use the restart data to convert the cold start to a warm start, which is required for JEDI coupling
cp $INPUTDATA_ROOT/RESTART/* INPUT
cd INPUT
for file in 20210323.060000.*; do new=$(echo $file | cut -c 17-) && mv -v -- "$file" "$new" ; done
${SED} -i 's/3    22/3    23/g' coupler.res
cd ..

# update model_configure file to turn off write component and do a warm start for our regression test
${SED} -i 's/quilting:                .true./quilting:                .false./g' model_configure
${SED} -i 's/start_day:               22/start_day:               23/g' model_configure
${SED} -i '/quilting_restart/d' model_configure

#change ufs.configure to use only 6 cores (no write component)
${SED} -i 's/0 7/0 5/g' ufs.configure
${SED} -i '/ATM_omp/d' ufs.configure

#switch to warm_start and turn off checksum comparisons
${SED} -i 's/warm_start = .false./warm_start = .true./g' input.nml
${SED} -i 's/make_nh = .true./make_nh = .false./g' input.nml
${SED} -i 's/na_init = 1/na_init = 0/g' input.nml
${SED} -i 's/external_ic = .true./external_ic = .false./g' input.nml
${SED} -i 's/nggps_ic = .true./nggps_ic = .false./g' input.nml
${SED} -i 's/mountain = .false./mountain = .true./g' input.nml
${SED} -i '/atmos_model_nml/a\  ignore_rst_cksum = .true.' input.nml
${SED} -i '/fv_core_nml/a\  ignore_rst_cksum = .true.' input.nml
${SED} -i '/&fms_nml/i \
&fms_io_nml\
  checksum_required = .false.\
/' input.nml




