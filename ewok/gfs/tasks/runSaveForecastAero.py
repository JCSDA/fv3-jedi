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

# Date
fcdate = conf['fc']['date']
base = conf['experiment']['expid'] + '.fc.' + yamltools.jedifnformat(fcdate) + "."

# Loop over steps
for sstep in conf['fc']['fcout']:

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
        source_file=f'{base}{sstep}.coupler.res',
        file_type='coupler.res',
        member=member,
    )

    file_types = ['fv_tracer.res']
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
                source_file=f'{base}{sstep}.{file_type}.tile{tile}.nc',
                file_extension='nc',
                file_type=file_type,
                tile=tile,
                member=member,
            )