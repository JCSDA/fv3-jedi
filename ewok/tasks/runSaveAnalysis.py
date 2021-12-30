#! /usr/bin/env python3

# (C) Copyright 2020-2021 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import sys
import os
import yamltools
import r2d2

conf = yamltools.configure_runtime(sys.argv[1])

# Check for working directory
if not os.path.exists(conf['currentdir']):
    raise RuntimeError('Working directory does not exist')
os.chdir(conf['currentdir'])

# Date
andate = conf['date']
base = conf['experiment']['expid'] + '.an.' + andate
filename = base + '.$(file_type).tile$(tile).nc'
cplrfile = base + '.coupler.res'

print("saveAnalysisRun filename = ", filename)

r2d2.store(
    model=conf['experiment']['model'],
    type='an',
    experiment=conf['experiment']['expid'],
    resolution=conf['resolution'],
    date=andate,
    source_file=filename,
    file_format='netcdf',
    file_type=['fv_core.res', 'fv_srf_wnd.res', 'fv_tracer.res', 'sfc_data'],
    tile=[1, 2, 3, 4, 5, 6],
)

print("saveAnalysisRun cplrfile = ", cplrfile)

r2d2.store(
    model='gfs_metadata',
    type='an',
    experiment=conf['experiment']['expid'],
    resolution=conf['resolution'],
    date=andate,
    source_file=cplrfile,
    file_type=['coupler.res'],
)
