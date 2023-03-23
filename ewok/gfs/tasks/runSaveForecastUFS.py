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

# Loop over steps, only saving the ones written out by UFS (contain RESTART in their name)
for sstep in conf['fc']['fcout']:
    if 'RESTART' in sstep:
        sstep = sstep.replace(" RESTART", "")
        final_date = yamltools.parse_datetime(fcdate) + yamltools.parse_timedelta(sstep)

        base = 'RESTART/' + final_date.strftime("%Y%m%d.%H%M%S")

        member = R2D2Data.DEFAULT_INT_VALUE
        if 'member' in conf:
            member = conf['member']

        R2D2Data.store(
            model=conf['experiment']['model'],
            item='forecast',
            step=sstep,
            experiment=conf['experiment']['expid'],
            resolution=conf['resolution'],
            date=fcdate,
            file_extension='',
            source_file=f'{base}.coupler.res',
            file_type='coupler.res',
            member=member,
        )

        file_types = ['fv_core.res', 'fv_srf_wnd.res', 'fv_tracer.res', 'sfc_data', 'phy_data']
        tiles = [1, 2, 3, 4, 5, 6]

        for file_type in file_types:
            for tile in tiles:
                R2D2Data.store(
                    model=conf['experiment']['model'],
                    item='forecast',
                    experiment=conf['experiment']['expid'],
                    step=sstep,
                    resolution=conf['resolution'],
                    date=fcdate,
                    source_file=f'{base}.{file_type}.tile{tile}.nc',
                    file_extension='nc',
                    file_type=file_type,
                    tile=tile,
                    member=member,
                )
