#! /usr/bin/env python3

# (C) Copyright 2020-2021 UCAR
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

# Define experiment to read from, current experiment by default
exp_read = conf['experiment']['expid']
if 'exp_source' in conf:
    exp_read = conf['exp_source']

# Fetch state

base = conf['experiment']['expid'] + '.fc.'
sdate = yamltools.jedifnformat(conf['fcdate']) + '.' + conf['fcstep']

fcstep = yamltools.parse_timedelta(conf['fcstep'])
if 'hack_step_bg' in conf and conf['hack_step_bg'] == True:
    pt3h = yamltools.parse_timedelta('PT3H')
    fcstep -= pt3h
hackstep = yamltools.jediformat(fcstep)

member = R2D2Data.DEFAULT_INT_VALUE
if 'member' in conf:
    member = conf['member']

R2D2Data.fetch(
    model=conf['experiment']['model'],
    item='forecast',
    step=hackstep,
    experiment=exp_read,
    resolution=conf['resolution'],
    date=conf['fcdate'],
    target_file=f'{base}{sdate}.coupler.res',
    file_extension='',
    file_type='coupler.res',
    member=member,
)

file_types = ['fv_core.res', 'fv_srf_wnd.res', 'fv_tracer.res', 'sfc_data']
tiles = [1, 2, 3, 4, 5, 6]

for file_type in file_types:
    for tile in tiles:
        R2D2Data.fetch(
            model=conf['experiment']['model'],
            item='forecast',
            experiment=exp_read,
            step=hackstep,
            resolution=conf['resolution'],
            date=conf['fcdate'],
            target_file=f'{base}{sdate}.{file_type}.tile{tile}.nc',
            file_extension='nc',
            file_type=file_type,
            tile=tile,
            member=member,
        )
