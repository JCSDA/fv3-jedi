#! /usr/bin/env python3

# (C) Copyright 2020-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import sys
import os
import yamltools
from r2d2 import fetch

conf = yamltools.configure_runtime(sys.argv[1])

# Check for working directory
if not os.path.exists(conf['currentdir']):
    raise RuntimeError('Working directory does not exist')
os.chdir(conf['currentdir'])

# Define experiment to read from, current experiment by default
exp_read = conf['experiment']['expid']
if 'exp_source' in conf:
    exp_read = conf['exp_source']

# Date and step
date = yamltools.parse_datetime(conf['date'])
offset = yamltools.parse_timedelta(conf['offset'])
fcdate = date - offset
fcstep = yamltools.parse_timedelta(conf['fcstep'])
vdate = fcdate + fcstep

# Fetch state

base = conf['experiment']['expid'] + '.fc.'
sdate = yamltools.jediformat(fcdate) + '.' + yamltools.jediformat(fcstep)
filename = base + sdate + '.$(file_type).tile$(tile).nc'
cplrfile = base + sdate + '.coupler.res'

fetch(
    model=conf['experiment']['model'],
    type='fc',
    experiment=exp_read,
    resolution=conf['resolution'],
    date=yamltools.jediformat(fcdate),
    step=conf['fcstep'],
    target_file=filename,
    file_format='netcdf',
    file_type=['fv_core.res', 'fv_srf_wnd.res', 'fv_tracer.res', 'sfc_data'],
    tile=[1, 2, 3, 4, 5, 6],
    fc_date_rendering='analysis',
)

fetch(
    model='gfs_metadata',
    type='fc',
    experiment=exp_read,
    resolution=conf['resolution'],
    date=yamltools.jediformat(fcdate),
    step=conf['fcstep'],
    target_file=cplrfile,
    file_type=['coupler.res'],
    fc_date_rendering='analysis',
)
