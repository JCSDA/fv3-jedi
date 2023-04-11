#! /usr/bin/env python3

# (C) Copyright 2020-2022 UCAR
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
andate = conf['an']['datetime']
base = conf['experiment']['expid'] + '.an.' + yamltools.jedifnformat(andate)

R2D2Data.store(
    model=conf['experiment']['model'],
    item='analysis',
    experiment=conf['experiment']['expid'],
    resolution=conf['resolution'],
    date=andate,
    source_file=f'{base}.bkg.nc',
    file_extension='nc',
    file_type='bkg',
)
