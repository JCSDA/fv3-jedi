#! /usr/bin/env python3

# (C) Copyright 2020-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import sys
import os
import yamltools
import r2d2

conf = yamltools.configure_runtime(sys.argv[1])

# Check for working directory
if not os.path.exists(conf['workdir']):
    raise RuntimeError('Working directory does not exist')
os.chdir(conf['workdir'])

# Date
fcdate = conf['fc']['date']
base = conf['experiment']['expid'] + '.fc.'

# Loop over steps
for sstep in conf['fc']['fcout']:
    print("saveForecastRun step = ", sstep)

    fcstep = yamltools.parse_timedelta(sstep)
    forecast_date = yamltools.parse_datetime(fcdate) + fcstep
    file_date = yamltools.jediformat(forecast_date)

    print("saveForecastRun for time = ", file_date)

    filename = base + file_date + '.$(file_type).nc'

    print("saveForecastRun filename = ", filename)

    r2d2.store(
        model=conf['experiment']['model'],
        type='fc',
        experiment=conf['experiment']['expid'],
        resolution=conf['resolution'],
        date=fcdate,
        step=sstep,
        source_file=filename,
        file_format='netcdf',
        file_type=['bkg'],
        fc_date_rendering='analysis',
    )
