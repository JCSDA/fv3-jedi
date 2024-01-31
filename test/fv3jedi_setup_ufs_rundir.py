#!/usr/bin/env python3
#
# Usage:
# fv3jedi_setup_ufs_rundir.py <UFS rundir> <test rundir>
#
# Because this script is lightweight (symlinks and mkdir), there is no
# check of whether the setup was previously done and could be skipped.

import os
import sys


ufs_rundir = sys.argv[1]
test_rundir = sys.argv[2]

# sanity checks
if not os.path.exists(ufs_rundir):
    raise Exception("Could not find UFS model rundir at path: " + ufs_rundir)
if not os.path.exists(test_rundir):
    raise Exception("Could not find dir for UFS-JEDI tests at path: " + test_rundir)

# symlink files from UFS rundir into test rundir
for item in os.listdir(ufs_rundir):
    # use absolute path to avoid symlink relative path shenanigans
    abs_ufs_path = os.path.abspath(os.path.join(ufs_rundir, item))
    test_path = os.path.join(test_rundir, item)
    if not os.path.exists(test_path):
        os.symlink(abs_ufs_path, test_path)
#   else:
#       # file exists; this is ok if the symlink was previously created, but otherwise it's an error
#       if not (os.path.islink(test_path) and (os.readlink(test_path) == abs_ufs_path)):
#           raise Exception("Could not link to UFS model rundir, path already exists: " + test_path)

# make new RESTART dir within rundir
new_restart_dir = os.path.join(test_rundir, 'RESTART')
if not os.path.exists(new_restart_dir):
    os.mkdir(new_restart_dir)
