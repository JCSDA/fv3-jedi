#!/usr/bin/env python3

#
# (C) Copyright 2020 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#

import click
import glob
from netCDF4 import Dataset
import numpy as np

# --------------------------------------------------------------------------------------------------
## @package rename_cold_starts
#  This function appends _cold to the name of the variables in the cold starts
#
#  e.g. t -> t_cold
#
# --------------------------------------------------------------------------------------------------

def abort(message):
    print('\n ABORT: '+message+'\n')
    raise(SystemExit)

# --------------------------------------------------------------------------------------------------

@click.command(help='Example: python rename_cold_starts.py 20200101.000000.gfs_data.tile*nc',
               context_settings={"ignore_unknown_options": True})
@click.argument('files', nargs=-1, type=click.Path())
@click.option('--user_prompt', help='Prompt before overwriting [true]', default=True, type=bool)
def main(files, user_prompt):

  # Print
  print('\n Running rename_cold_starts on the following files:')
  for file in files:
    print('  '+file)
  print(' ')

  # Prompt user to avoid accidental overwrite
  if user_prompt:
    yesno = input("This will overwrite the above files, continue? [yes]/no ") or "yes"
    if not yesno == 'yes':
      abort("User did not select yes to continue option.")

  # Fields to append with _cold
  fields = ['ps','w','zh','t','delp','sphum','liq_wat','o3mr','ice_wat','rainwat','snowwat',
            'graupel','u_w','v_w','u_s','v_s']

  # Rename fields in the netCDF files
  # ---------------------------------
  for file in files:

    ncfile = Dataset(file, mode='r+')

    for field in fields:

      try:
        ncfile.renameVariable(field,field+'_cold')
      except:
        print('Field '+field+' could not be renamed in '+file)

    ncfile.close()

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    main()

# --------------------------------------------------------------------------------------------------
