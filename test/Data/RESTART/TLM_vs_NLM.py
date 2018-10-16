from netCDF4 import Dataset
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fcst_nl_ana = 


#User input required for the follwing:
plot_diff = 0         #Plot path1/file - path2/file
model = 'gfs'
cube = 48
filetype = 'png'

path1  = './'                         #Path of first/only file
file_tplt_befr1 = '20181004.000000.Jnorm.fv_tracer.res.tile'  #Filename befor tile number
file_tplt_aftr = '.nc'                #Filename after tile number

if (cube == 48):
    if (model == 'geos'):
        path2  = '../INPUTS/GEOS_c48/'
        file_tplt_befr2 = '20180415.000000.fv_core.res.tile'
    elif (model == 'gfs'):
        path2  = '../INPUTS/GFS_c48/ENSEMBLE/mem001/RESTART/'
        file_tplt_befr2 = '20180415.000000.fv_core.res.tile'
elif (cube == 96):
    path2  = ''
    file_tplt_befr2 = 'fv_core.res.tile'

xdimvar = 'xaxis_1'                  #What to read to get dimension
ydimvar = 'yaxis_1'                  #What to read to get dimension
zdimvar = 'zaxis_1'                  #What to read to get dimension
readvar = 'sphum'                    #Variable to plot
Dim2dor3d = '3D'                     #Is this 2D or 3D field?
plot_level = 50                      #If 3D plot this level

readvarlat = 'grid_lat'              #Variable to plot
readvarlon = 'grid_lon'              #Variable to plot

#Tile 1 handle for getting dimension
fh1 = Dataset(path1 + file_tplt_befr1 + str(1) + file_tplt_aftr, mode='r')
fh2 = Dataset(path2 + file_tplt_befr2 + str(1) + file_tplt_aftr, mode='r')

#Dimensions
npx = len(fh1.dimensions[xdimvar])
npy = len(fh1.dimensions[ydimvar])
npz = 0

npxr = npx
npyr = npy

if (readvar == 'u'):
   npyr = npyr + 1
if (readvar == 'v'):
   npxr = npxr + 1

if Dim2dor3d == '3D':
   npz = len(fh1.dimensions[zdimvar])
   fr = np.zeros(6*npz*npyr*npxr).reshape(6,npz,npyr,npxr)
else:
   fr = np.zeros(6*npyr*npxr).reshape(6,npyr,npxr)

flat = np.zeros(6*npy*npx).reshape(6,npy,npx)
flon = np.zeros(6*npy*npx).reshape(6,npy,npx)

#Read file
for tile in range(6):
    file_tile1 = file_tplt_befr1 + str(tile+1) + file_tplt_aftr
    file_tile2 = file_tplt_befr2 + str(tile+1) + file_tplt_aftr
    pathfile1 = path1 + file_tile1
    pathfile2 = path2 + file_tile2
    print(pathfile1)
    if plot_diff == 1:
        print(' '+pathfile2)
    fh1 = Dataset(pathfile1, mode='r')
    fh2 = Dataset(pathfile2, mode='r')
    if plot_diff == 1:
        fr[tile,:,:,:] = fh1.variables[readvar][:] - fh2.variables[readvar][:]
    else:
        fr[tile,:,:,:] = fh1.variables[readvar][:]
