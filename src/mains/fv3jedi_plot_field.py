#!/usr/bin/env python3
#
# (C) Copyright 2020 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import cartopy.crs as ccrs
import click
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import os

# --------------------------------------------------------------------------------------------------
## @package fv3jedi_plot_field.py
#  This is a utility for plotting an fv3jedi field that has been output on lon/lat grid
#
#  1. Run the fv3jedi_convertstate.x application outputting a lonlat set of fields. See
#     fv3-jedi/test/testinput/convertstate_gfs_c2ll.yaml for an example of doing this.
#
#  2. Run this application with the inputs corresponding to the field and level that you
#     wish to plot.
#
# --------------------------------------------------------------------------------------------------

def abort(message):
    print('\n ABORT: '+message+'\n')
    raise(SystemExit)

# --------------------------------------------------------------------------------------------------

@click.command()
@click.option('--inputfile', required=True,  help='NetCDF input containing Lon/Lat fields')
@click.option('--fieldname', required=True,  help='Name of field to plot')
@click.option('--layer',     required=False, help='Model layer to plot. [reqiured if 3D variable]', type=int)
@click.option('--showfig',   required=False, help='Display figure if set to true', default=False)
def main(inputfile, fieldname, layer, showfig):

    # Open the file
    print('\nOpening ', inputfile, 'for reading')
    ncfile = netCDF4.Dataset(inputfile, mode='r')

    # Get metadata from the file
    npx = ncfile.dimensions["lon"].size
    npy = ncfile.dimensions["lat"].size
    npz = ncfile.dimensions["lev"].size
    lons = ncfile.variables["lons"][:]
    lats = ncfile.variables["lats"][:]

    # Print field dimensions
    print(" Grid dimensions", npx, 'x', npy, 'x', npz)

    # Get field units from the file
    units = ncfile.variables[fieldname].units

    # Zero out array to fill with field
    field = np.zeros((npy, npx))

    # Check if field is two or three dimensions
    if len(ncfile.variables[fieldname].shape) == 4:

      # User must provide layer/level to plot if 3D
      if (layer == None):
          abort("If plotting 3D variable user must provide layer with --layer")

      # Message and read the field at provided layer
      print(" Reading layer ", layer, " from field ", fieldname)
      field[:,:] = ncfile.variables[fieldname][:,layer-1,:,:]

      # Set plot title and output file to include level plotted
      title = "Contour of "+fieldname+" ("+units+") for layer "+str(layer)
      outfile = os.path.splitext(inputfile)[0]+"_"+fieldname+"_layer-"+str(layer)+".png"

    elif len(ncfile.variables[fieldname].shape) == 3:

      # Message and read the field at provided layer
      print(" Reading field ", fieldname)
      field[:,:] = ncfile.variables[fieldname][:,:]
      title = "Contour of "+fieldname+" ("+units+")"
      outfile = os.path.splitext(inputfile)[0]+"_"+fieldname+".png"

    # Close the file
    ncfile.close()

    # Check if field has positve and negative values
    # ----------------------------------------------
    if np.min(field) < 0:
      cmax = np.max(np.abs(field))
      cmin = -cmax
      cmap = 'RdBu'
    else:
      cmax = np.max(field)
      cmin = np.min(field)
      cmap = 'nipy_spectral'

    levels = np.linspace(cmin,cmax,25)

    # Create two dimensional contour plot of field
    # --------------------------------------------

    # Set the projection
    projection = ccrs.PlateCarree()

    # Create figure to hold plot
    fig = plt.figure(figsize=(10, 5))

    # Just one subplot for now
    ax = fig.add_subplot(1, 1, 1, projection=projection)

    # Contour the field
    im = ax.contourf(lons, lats, field,
                     transform=projection,
                     cmap=cmap,
                     levels=levels)

    # Add coast lines to the plot
    ax.coastlines()

    # Add labels to the plot
    ax.set_xticks(np.linspace(-180, 180, 5), crs=projection)
    ax.set_yticks(np.linspace(-90, 90, 5), crs=projection)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitide')
    ax.set_title(title)
    ax.set_global()

    # Add a colorbar for the filled contour.
    fig.colorbar(im)

    # Show the figure
    print(" Saving figure as", outfile, "\n")
    plt.savefig(outfile)
    if (showfig):
        plt.show()


# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    main()

# --------------------------------------------------------------------------------------------------

