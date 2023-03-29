#! /usr/bin/env python3

# (C) Copyright 2023 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import sys
import os
import yamltools
import shutil

conf = yamltools.parse_config(sys.argv[1])

# Check for working directory
if not os.path.exists(conf['workdir']):
    raise RuntimeError('Working directory does not exist')
os.chdir(conf['workdir'])

# Harvest info from config:
fv3repo = conf['fv3repo']
static_data = conf['static_data']
fixdir = conf['fixdir']
resol = conf['resol']
nlevs = str(conf['nlevs'])
layout = conf['layout']

# Files that need to be at the top-level of ewok/<expid>/
top_src = [
fv3repo+'/test/Data/fv3files/akbk'+nlevs+'.nc4',
fv3repo+'/test/Data/fieldmetadata/ufs.yaml',
fv3repo+'/test/Data/fv3files/field_table_ufs'
]

top_dest = [
fixdir+'/akbk'+nlevs+'.nc4',
fixdir+'/ufs.yaml',
fixdir+'/field_table'
]

for ii in range(len(top_src)):
    if not os.path.isfile(top_dest[ii]):
        shutil.copy(top_src[ii], top_dest[ii])

# Copy all fix files needed to run the model in workdir/<expid>/<cycle>/UFSFixFiles
# the path to this repo needs to be in `model.ufs_run_directory`

# Make UFSFiles dir within workdir
ufsfiles_dir = os.path.join(fixdir, 'UFSFixFiles')
if not os.path.exists(ufsfiles_dir):
    os.mkdir(ufsfiles_dir)

release = conf['listfiles']['release']
cycle_src = conf['listfiles']['fix files']

for ii in range(len(cycle_src)):
    src = static_data+'/'+release+'/ufs/'+resol+'/'+cycle_src[ii]
    dest = ufsfiles_dir+'/'+cycle_src[ii]
    if not os.path.isfile(dest):
        shutil.copy(src, dest)

# Make UFSFiles/INPUT dir within workdir
ufsfiles_input_dir = os.path.join(ufsfiles_dir, 'INPUT')
if not os.path.exists(ufsfiles_input_dir):
    os.mkdir(ufsfiles_input_dir)

cycle_input_src = conf['listfiles']['input fix files']

for ii in range(len(cycle_input_src)):
    src = static_data+'/'+release+'/ufs/'+resol+'/INPUT/'+cycle_input_src[ii]
    dest = ufsfiles_input_dir+'/'+cycle_input_src[ii]
    if not os.path.isfile(dest):
        shutil.copy(src, dest)

# Resolution dependent files:
res_files_src = [
fv3repo+'/test/Data/fv3files/input_ufs_'+resol+'_'+layout+'.nml',
fv3repo+'/test/Data/fv3files/input_ufs_'+resol+'_'+layout+'.nml',
fv3repo+'/test/Data/fv3files/nems.configure.'+layout
]

res_files_dest = [
fixdir+'/input.nml',
fixdir+'/UFSFixFiles/input.nml',
fixdir+'/UFSFixFiles/nems.configure'
]

for ii in range(len(res_files_src)):
    if not os.path.isfile(res_files_dest[ii]):
        shutil.copy(res_files_src[ii], res_files_dest[ii])
