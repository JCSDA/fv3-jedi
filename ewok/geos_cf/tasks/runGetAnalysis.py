#! /usr/bin/env python3

# (C) Copyright 2020-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import sys
import os
import yamltools
from r2d2 import fetch

conf = yamltools.configure_runtime(sys.argv[1])

# Check for working directory
if not os.path.exists(conf['workdir']):
    raise RuntimeError('Working directory does not exist')
os.chdir(conf['workdir'])

# Define experiment to read from, current experiment by default
exp_read = conf['experiment']['expid']
if 'exp_source' in conf:
    exp_read = conf['exp_source']

# Date
andate = yamltools.parse_datetime(conf['date'])

# Fetch state

base = conf['experiment']['expid'] + '.an.'
sdate = yamltools.jediformat(andate)
filename = base + sdate + '.$(file_type).nc'

fetch(
    model=conf['experiment']['model'],
    type='an',
    experiment=exp_read,
    resolution=conf['resolution'],
    date=yamltools.jediformat(andate),
    target_file=filename,
    file_format='netcdf',
    file_type=['bkg'],
    fc_date_rendering='analysis',
)
