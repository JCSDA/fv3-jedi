#!/usr/bin/env python3
"""
@author: Benjamin Menetrier
@description: plotting facility for FV3
"""

# (C) Copyright 2021-2022 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

import cartopy.crs as ccrs
import copy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from netCDF4 import Dataset
import numpy as np
import os
import time
import subprocess
import sys
import yamltools

# -----------------------------------------------------------------------------

# Initial time
initial_time = time.perf_counter()

# -----------------------------------------------------------------------------
# Harvest values from the config

conf = yamltools.configure_runtime(sys.argv[1])

filepath = conf['filepath']
variables = conf['variables']
levels = conf['levels']
output = conf['output']
gridfiledir = conf['gridfiledir']

# basefilepath ?
# -----------------------------------------------------------------------------

# Lon/lat of view
lonview = -45.0
latview = 45.0

# -----------------------------------------------------------------------------
# Set default string values

if 'colormap' in conf:
    colormap = conf['colormap']
else:
    if 'centered' in conf:
        colormap = "seismic"
    else:
        colormap = "viridis"

centered = True if 'centered' in conf else False
diagnostic = conf['diagnostic'] if 'diagnostic' in conf else 'default'


# -----------------------------------------------------------------------------
# Print arguments

print("Parameters:")
print(" - gridfiledir: " + gridfiledir)
print(" - filepath: " + filepath)
print(" - variable:")
print(variables)
print(" - levels:")
print(levels)
print(" - diagnostic: " + diagnostic)
print(" - colormap: " + colormap)
print(" - centered: " + str(centered))
print(" - output: " + output)

# -----------------------------------------------------------------------------
# Convert levels list to array of integers
levels = np.array(list(map(int, levels)))
# -----------------------------------------------------------------------------

# Check file extension
if not filepath.endswith(".nc"):
    print("   Error: filepath extension should be .nc")
    sys.exit(1)

nvars = len(variables)

units = [None] * nvars
long_name = [None] * nvars
nx = [None] * nvars
ny = [None] * nvars
fld = [None] * nvars
vmin = [None] * nvars
vmax = [None] * nvars

for itile in range(0, 6):
    # Open data file
    filename = filepath.replace(".nc", ".tile" + str(itile+1) + ".nc")
    fdata = Dataset(filename, "r", format="NETCDF4")

    # Read field
    for var in range(nvars):

        fld_tmp = fdata[variables[var]][0,levels-1,:,:]
        units[var] = fdata[variables[var]].units
        long_name[var] = fdata[variables[var]].long_name

        if itile == 0:
            # Get shape
            shp = np.shape(fld_tmp)
            nz = shp[0]
            ny[var] = shp[1]
            nx[var] = shp[2]
            grid = min(ny[var], nx[var])

            # Initialize field
            fld[var] = np.zeros((6, nz, ny[var], nx[var]))

        # Copy field
        fld[var][itile,:,:,:] = fld_tmp

        # Compute min/max
        if centered:
            vmax[var] = np.max(np.abs(fld[var]))
            vmin[var] = -vmax
        else:
            vmin[var] = np.min(fld[var])
            vmax[var] = np.max(fld[var])


# Open grid file and read lons/lats + levels to pressure
fgrid = Dataset(gridfiledir + "/fv3grid_c" + str(grid).zfill(4) + ".nc4",
                "r", format="NETCDF4")
vlons = np.degrees(fgrid["vlons"][:,:,:])
vlats = np.degrees(fgrid["vlats"][:,:,:])

# Compute all pressures in hPa: P = ak + bk Ps
# I got a point in the water for the Ps and we are looking at levels 850hPa and less
akbk = Dataset(gridfiledir + "/akbk127.nc4", "r", format="NETCDF4")
ak = akbk["ak"]
bk = akbk["bk"]
pressure = [None] * np.shape(ak)[0]
for ll in range(0, np.shape(ak)[0]):
    pressure[ll] = (ak[ll] + bk[ll] * 101300)/100

for iz in range(0, nz):

    # Projection
    projection = ccrs.Orthographic(lonview, latview)

    for var in range(nvars):
        # Figure title
        title = long_name[var] + " (" + units[var] + ") at pressure " + str(round(pressure[levels[iz]], 0)) + " hPa - C" + str(nx[var])
        title = title.replace("_", " ")

        # Output
        outputname = output + "_" + variables[var] + "_" + str(levels[iz])

        # Initialize figure
        fig,ax = plt.subplots(figsize=(8,8),subplot_kw=dict(projection=projection))
        ax.set_global()
        ax.coastlines()

        # Colormap
        cmap = copy.copy(cm.get_cmap(colormap))
        cmap.set_bad('gray', 1)

        # Normalization
        norm = plt.Normalize(vmin=vmin[var], vmax=vmax[var])

        # Figure title
        plt.title(title)
        nx[var] = min(nx[var], ny[var])
        ny[var] = nx[var]

        # Loop over tiles
        for itile in range(0, 6):
            # Loop over polygons
            for iy in range(0, ny[var]):
                for ix in range(0, nx[var]):
                    # Polygon coordinates
                    xy = [[vlons[itile,iy+0,ix+0], vlats[itile,iy+0,ix+0]],
                          [vlons[itile,iy+0,ix+1], vlats[itile,iy+0,ix+1]],
                          [vlons[itile,iy+1,ix+1], vlats[itile,iy+1,ix+1]],
                          [vlons[itile,iy+1,ix+0], vlats[itile,iy+1,ix+0]]]

                    # Add polygon
                    ax.add_patch(mpatches.Polygon(xy=xy, closed=True, facecolor=cmap(norm(fld[var][itile,iz,iy,ix])),transform=ccrs.Geodetic()))

        # Set colorbar
        sm = cm.ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])
        plt.colorbar(sm, orientation="vertical",shrink=0.8)

        # Save and close figure
        plt.savefig(outputname + ".png", format="png", dpi=300)
        plt.close()
        print("Created ", outputname+".png")

        # Trim figure with mogrify if available
        info = subprocess.getstatusoutput('mogrify -help')
        if info[0] == 0:
            subprocess.run(["mogrify", "-trim", output + ".png"])

# -----------------------------------------------------------------------------

# Final time
final_time = time.perf_counter()

# Print timing
print(f"plot_gfs.py executed in {final_time - initial_time:0.4f} seconds")
