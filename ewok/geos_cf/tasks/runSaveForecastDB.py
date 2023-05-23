#! /usr/bin/env python3

# (C) Copyright 2020-2023 UCAR
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

# Date
fcdate = conf['fc']['date']

# Filename
expid = conf['experiment']['expid']
base = expid + '.' + conf['resolution']+ '.fc.'

# Target name
target_expid = conf['target_expid']

# Loop over steps
for sstep in conf['fc']['fcout']:

    fcstep = yamltools.parse_timedelta(sstep)
    forecast_date = yamltools.parse_datetime(fcdate) + fcstep
    file_date = yamltools.jedifnformat(forecast_date)

    if 'target_DB' in conf:
        R2D2Data.store(
            model=conf['experiment']['model'],
            item='forecast',
            step=sstep,
            experiment=conf['experiment']['expid'],
            resolution=conf['resolution'],
            date=fcdate,
            source_file=f'{base}{file_date}.bkg.nc',
            file_extension='nc',
            file_type='bkg',
            data_store=conf['target_DB'],
        )

    else:
        R2D2Data.store(
            model=conf['experiment']['model'],
            item='forecast',
            step=sstep,
            experiment=conf['experiment']['expid'],
            resolution=conf['resolution'],
            date=fcdate,
            source_file=f'{base}{file_date}.bkg.nc',
            file_extension='nc',
            file_type='bkg',
        ) 