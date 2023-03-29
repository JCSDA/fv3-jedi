#! /usr/bin/env python3

# (C) Copyright 2023 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import sys
import os
import yamltools
from r2d2 import R2D2Data

conf = yamltools.configure_runtime(sys.argv[1])

# Check for working directory
if not os.path.exists(conf['workdir']):
    raise RuntimeError('Working directory does not exist')
os.chdir(conf['workdir'])

# Fetch input state from r2d2 database
# ------------------------------------
date = conf['andate']
exp_read = conf['experiment']['expid']
if 'exp_source' in conf:
    exp_read = conf['exp_source']

base = 'INPUT/'
basedir = os.path.join(conf['workdir'], base)
if not os.path.exists(basedir):
    os.mkdir(basedir)

member = R2D2Data.DEFAULT_INT_VALUE
if 'member' in conf:
    member = conf['member']

R2D2Data.fetch(
    model=conf['experiment']['model'],
    item='analysis',
    experiment=exp_read,
    resolution=conf['resolution'],
    date=date,
    target_file=f'{base}/coupler.res',
    file_type='coupler.res',
    file_extension='',
    member=member,
)

file_types = ['fv_core.res', 'fv_srf_wnd.res', 'fv_tracer.res', 'sfc_data', 'phy_data']
tiles = [1, 2, 3, 4, 5, 6]

for file_type in file_types:
    for tile in tiles:
        R2D2Data.fetch(
            model=conf['experiment']['model'],
            item='analysis',
            experiment=exp_read,
            resolution=conf['resolution'],
            date=date,
            target_file=f'{base}/{file_type}.tile{tile}.nc',
            file_extension='nc',
            file_type=file_type,
            tile=tile,
            member=member,
        )


# Create links from static files dir to current workdir
# -----------------------------------------------------

fcworkdir = conf['fcworkdir']
ufs_files = conf['ufs_modeldir']

print("Debug: ufs_files read from:", ufs_files)
print("Debug: fcworkdir current workdir", fcworkdir)

list_ufs_files = os.listdir(ufs_files)

for item in list_ufs_files:
    abs_ufs_path = os.path.join(ufs_files, item)
    workdir_path = os.path.join(fcworkdir, item)

    if os.path.isfile(abs_ufs_path) and not os.path.exists(workdir_path):
        os.symlink(abs_ufs_path, workdir_path)

    elif os.path.isdir(abs_ufs_path):
        if not os.path.isdir(workdir_path):
            os.mkdir(workdir_path)
        for input_item in os.listdir(abs_ufs_path):
            input_abs_ufs_path = os.path.join(abs_ufs_path, input_item)
            input_workdir_path = os.path.join(workdir_path, input_item)
            if not os.path.exists(input_workdir_path):
                os.symlink(input_abs_ufs_path, input_workdir_path)

# Make new RESTART dir within workdir
new_restart_dir = os.path.join(fcworkdir, 'RESTART')
if not os.path.exists(new_restart_dir):
    os.mkdir(new_restart_dir)


# Set up model_configure
# ----------------------
fv3repo = conf['fv3repo']+'/test/Data/fv3files'
modelconfig_read_path = os.path.join(fv3repo, 'model_configure_'+conf['resolution'])
modelconfig_write_path = os.path.join(fcworkdir, 'model_configure')

# Harvest forecast length and output frequency
fc_length = yamltools.parse_timedelta(conf['fc_length'])
fc_length = str( fc_length.days*24 + fc_length.seconds//3600 )

fc_freq = yamltools.parse_timedelta(conf['fc_freq'])
fc_freq = str( fc_freq.days*24 + fc_freq.seconds//3600 )

# For now we read info from coupler.res
couplerfile = os.path.join(fcworkdir, 'INPUT/coupler.res')
coupler = open(couplerfile, 'r')

content = coupler.readlines()
values = content[1].split()
start_year = values[0]
start_month = values[1].zfill(2)
start_day = values[2].zfill(2)
start_hour = values[3].zfill(2)

values = content[0].split()
calendar_list = ['no_calendar', 'thirty_day_months', 'julian', 'gregorian', 'noleap']
calendar = calendar_list[int(values[0])]

coupler.close()

# Write model_config
model_config = open(modelconfig_write_path, 'w')

with open(modelconfig_read_path, 'r') as model_config_tmpl:
    for line in model_config_tmpl:
        line2 = line.replace('%year%', start_year)
        line2 = line2.replace('%month%', start_month)
        line2 = line2.replace('%day%', start_day)
        line2 = line2.replace('%hour%', start_hour)
        line2 = line2.replace('%fc_length%', fc_length)
        line2 = line2.replace('%fc_freq%', fc_freq)
        line2 = line2.replace('%calendar%', calendar)
        model_config.write(line2)
model_config.close()
