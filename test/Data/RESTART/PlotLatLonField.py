from netCDF4 import Dataset
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic
from matplotlib.ticker import FixedLocator
from matplotlib import rcParams

#User input required for the follwing:
plot_diff = 0         #Plot path1/file - path2/file
model = 'gfs'
cube = 180
filetype = 'png'

file_tplt_aftr = '.nc'                #Filename after tile number

datetime = '20180414_210000'

file2 = 'hyb-3DVar-fgat.rsonde.c48.latlon.'+datetime+'z.nc4'
file1 = 'forecast.fv3.c48.latlon.'+datetime+'z.nc4'

xdimvar = 'lon'                  #What to read to get dimension
ydimvar = 'lat'                  #What to read to get dimension
zdimvar = 'lev'                  #What to read to get dimension
readvar = 't'                    #Variable to plot
Dim2dor3d = '3D'                     #Is this 2D or 3D field?
plot_level = 60                      #If 3D plot this level

print(file1)
print(file2)

fh1 = Dataset(file1, mode='r')
fh2 = Dataset(file2, mode='r')

#Dimensions
im = len(fh1.dimensions[xdimvar])
jm = len(fh1.dimensions[ydimvar])
lm = len(fh1.dimensions[zdimvar])

lons = fh1.variables['lons'][:] - 180.0
lats = fh1.variables['lats'][:]

if Dim2dor3d == '3D':
   fr = np.zeros(lm*jm*im).reshape(lm,jm,im)
else:
   fr = np.zeros(jm*im).reshape(jm,im)

if plot_diff == 1:
    fr[:,:,:] = fh1.variables[readvar][:] - fh2.variables[readvar][:]
else:
    fr[:,:,:] = fh1.variables[readvar][:]


#Get level to plot
f = np.zeros(jm*im).reshape(jm,im)
if Dim2dor3d == '3D':
   f[:,:] = fr[plot_level,0:jm,0:im]
else:
   f = fr

#Contour levels
maxf = np.nanmax(np.abs(f))
minf = np.nanmin(f)

if minf < 0:
    minf = -maxf
incf = (maxf - minf)/101
clev = np.arange(minf,maxf+incf,incf)

incf = (maxf - minf)/10
ctic = np.arange(minf,maxf+incf,incf)

#Colormap
cmap = plt.cm.seismic

itit = '  '
if plot_diff == 1:
   itit = ' increment '

fig = plt.figure(figsize=(14,8))
ax = fig.add_axes([0.1,0.1,0.8,0.8])

m = Basemap(projection='cyl', \
            llcrnrlat=-90,    \
            urcrnrlat=90,     \
            llcrnrlon=-180,   \
            urcrnrlon=180,    \
            resolution='c',   \
            suppress_ticks   = False )

m.drawcoastlines()


fp, lonsp = addcyclic(f, lons) #Add overlap

lonsp, lats = np.meshgrid(lonsp, lats) #Mesh 
x, y = m(lonsp, lats) #Project to map

cp = m.contourf(x,y,fp,clev,cmap=plt.cm.seismic)

if Dim2dor3d == '3D':
   plt.title('Cubed sphere plot of ' + readvar + itit + 'at level: ' + str(plot_level))
else:
   plt.title('Cubed sphere plot of ' + readvar)

fig.patch.set_facecolor('grey')

plt.savefig('LatLonPlot_Field-'+readvar+'_Level-'+str(plot_level)+'.'+filetype, bbox_inches='tight',facecolor=fig.get_facecolor())
