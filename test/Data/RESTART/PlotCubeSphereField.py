from netCDF4 import Dataset
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap

#User input required for the follwing:
plot_diff = 1         #Plot path1/file - path2/file
model = 'fv3gfs'
cube = 48
filetype = 'png'

path1  = './'                         #Path of first/only file
file_tplt_befr1 = '20180415.000000.3D-Var.fv_core.res.tile'  #Filename befor tile number
file_tplt_aftr = '.nc'                #Filename after tile number

if (cube == 48):
    if (model == 'geos'):
        path2  = '../INPUTS/GEOS_c48/'
        file_tplt_befr2 = '20180415.000000.fv_core.res.tile'
    elif (model == 'fv3gfs'):
        path2  = '../INPUTS/FV3GFS_c48/ENSEMBLE/mem001/RESTART/'
        file_tplt_befr2 = '20180415.000000.fv_core.res.tile'
elif (cube == 96):
    path2  = ''
    file_tplt_befr2 = 'fv_core.res.tile'

xdimvar = 'xaxis_1'                  #What to read to get dimension
ydimvar = 'yaxis_2'                  #What to read to get dimension
zdimvar = 'zaxis_1'                  #What to read to get dimension
readvar = 'T'                        #Variable to plot
Dim2dor3d = '3D'                     #Is this 2D or 3D field?
plot_level = 40                      #If 3D plot this level

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
    flat[tile,:,:] = fh1.variables[readvarlat][:]
    flon[tile,:,:] = fh1.variables[readvarlon][:]
    if plot_diff == 1:
        fr[tile,:,:,:] = fh1.variables[readvar][:] - fh2.variables[readvar][:]
    else:
        fr[tile,:,:,:] = fh1.variables[readvar][:]


#Get level to plot
f = np.zeros(6*npy*npx).reshape(6,npy,npx)
if Dim2dor3d == '3D':
   f[:,:,:] = fr[:,plot_level,0:npy,0:npx]
else:
   f = fr

#Initial arrays
f12 = np.zeros(12*npy*npx).reshape(12,npy,npx)
fp = np.zeros(12*npy*npx).reshape(4*npy,3*npx)
f12lat = np.zeros(12*npy*npx).reshape(12,npy,npx)
f12lon = np.zeros(12*npy*npx).reshape(12,npy,npx)

#Transpose, rotate etc
f12lat[0,:,:] = np.rot90(flat[2,:,:],2)
f12lat[1,:,:] = np.rot90(flat[2,:,:],3)
f12lat[2,:,:] = flat[2,:,:]
f12lat[3,:,:] = np.rot90(flat[2,:,:],1)
f12lat[4,:,:] = np.fliplr(np.transpose(flat[0,:,:]))
f12lat[5,:,:] = np.fliplr(np.transpose(flat[1,:,:]))
f12lat[6,:,:] = flat[3,:,:]
f12lat[7,:,:] = flat[4,:,:]
f12lat[8,:,:] = np.rot90(flat[5,:,:],3)
f12lat[9,:,:] = np.rot90(flat[5,:,:],2)
f12lat[10,:,:] = np.rot90(flat[5,:,:],1)
f12lat[11,:,:] = flat[5,:,:]

f12lon[0,:,:] = np.rot90(flon[2,:,:],2)
f12lon[1,:,:] = np.rot90(flon[2,:,:],3)
f12lon[2,:,:] = flon[2,:,:]
f12lon[3,:,:] = np.rot90(flon[2,:,:],1)
f12lon[4,:,:] = np.fliplr(np.transpose(flon[0,:,:]))
f12lon[5,:,:] = np.fliplr(np.transpose(flon[1,:,:]))
f12lon[6,:,:] = flon[3,:,:]
f12lon[7,:,:] = flon[4,:,:]
f12lon[8,:,:] = np.rot90(flon[5,:,:],3)
f12lon[9,:,:] = np.rot90(flon[5,:,:],2)
f12lon[10,:,:] = np.rot90(flon[5,:,:],1)
f12lon[11,:,:] = flon[5,:,:]

f12[0,:,:] = np.rot90(f[2,:,:],2)
f12[1,:,:] = np.rot90(f[2,:,:],3)
f12[2,:,:] = f[2,:,:]
f12[3,:,:] = np.rot90(f[2,:,:],1)
f12[4,:,:] = np.fliplr(np.transpose(f[0,:,:]))
f12[5,:,:] = np.fliplr(np.transpose(f[1,:,:]))
f12[6,:,:] = f[3,:,:]
f12[7,:,:] = f[4,:,:]
f12[8,:,:] = np.rot90(f[5,:,:],3)
f12[9,:,:] = np.rot90(f[5,:,:],2)
f12[10,:,:] = np.rot90(f[5,:,:],1)
f12[11,:,:] = f[5,:,:]

#Cat to single array
fp[:,0*npx+0:npx] = np.concatenate([f12[0,:,:], f12[1,:,:], f12[2 ,:,:], f12[3 ,:,:]])
fp[:,1*npx:2*npx] = np.concatenate([f12[4,:,:], f12[5,:,:], f12[6 ,:,:], f12[7 ,:,:]])
fp[:,2*npx:3*npx] = np.concatenate([f12[8,:,:], f12[9,:,:], f12[10,:,:], f12[11,:,:]])
fp = np.flipud(np.transpose(fp))

#Contour levels
maxf = np.nanmax(np.abs(fp))
minf = np.nanmin(fp)

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


# #Plot with coast lines
# 
# fig = plt.figure(figsize=(14,8))
# 
# axN = fig.add_subplot(3,4,2)
# ax1 = fig.add_subplot(3,4,5)
# ax2 = fig.add_subplot(3,4,6)
# ax3 = fig.add_subplot(3,4,7)
# ax4 = fig.add_subplot(3,4,8)
# axS = fig.add_subplot(3,4,10)
# 
# bl = 48.0 #np.min(f12lat[0,:,:])
# 
# mapN = Basemap(projection='npstere', boundinglat= bl,lon_0 = np.min(0.5*(f12lon[5,npx/2-1,:] + f12lon[5,npx/2,:])), ax = axN)
# mapS = Basemap(projection='spstere', boundinglat=-bl,lon_0 = np.min(0.5*(f12lon[7,npx/2-1,:] + f12lon[7,npx/2,:])), ax = axS)
# 
# map1 = Basemap(projection='merc', llcrnrlat = np.min(f12lat[4,:,:]), urcrnrlat = np.max(f12lat[4,:,:]),\
#                                   llcrnrlon = -(360.0-np.min(f12lon[4,1,:])), urcrnrlon = np.max(f12lon[4,npx-1,:]), ax=ax1)
# map2 = Basemap(projection='merc', llcrnrlat = np.min(f12lat[5,:,:]), urcrnrlat = np.max(f12lat[5,:,:]),\
#                                   llcrnrlon = np.min(f12lon[5,:,:]), urcrnrlon = np.max(f12lon[5,:,:]), ax=ax2)
# map3 = Basemap(projection='merc', llcrnrlat = np.min(f12lat[6,:,:]), urcrnrlat = np.max(f12lat[6,:,:]),\
#                                   llcrnrlon = np.min(f12lon[6,:,:]), urcrnrlon = np.max(f12lon[6,:,:]), ax=ax3)
# map4 = Basemap(projection='merc', llcrnrlat = np.min(f12lat[7,:,:]), urcrnrlat = np.max(f12lat[7,:,:]),\
#                                   llcrnrlon = np.min(f12lon[7,:,:]), urcrnrlon = np.max(f12lon[7,:,:]), ax=ax4)
# 
# mapN.drawcoastlines()
# map1.drawcoastlines()
# map2.drawcoastlines()
# map3.drawcoastlines()
# map4.drawcoastlines()
# mapS.drawcoastlines()
# 
# x,y = np.meshgrid(np.linspace(mapN.xmin,mapN.xmax,npx),np.linspace(mapN.ymin,mapN.ymax,npx))
# csN = mapN.contourf(x,y,np.rot90(f12[1,:,:]),clev,cmap=cmap)
# 
# x,y = np.meshgrid(np.linspace(mapS.xmin,mapS.xmax,npx),np.linspace(mapS.ymin,mapS.ymax,npx))
# csS = mapS.contourf(x,y,np.rot90(f12[9,:,:],3),clev,cmap=cmap)
# 
# latplt = np.linspace(np.min(f12lat[4,:,:]), np.max(f12lat[4,:,:]), num=npx)
# lonplt = np.linspace(-(360.0-np.min(f12lon[4,1,:])), np.max(f12lon[4,npx-1,:]), num=npx)
# lonpltm, latpltm = np.meshgrid(lonplt, latplt)
# x, y = map1(lonpltm,latpltm)
# cs1 = map1.contourf(x,y,np.rot90(f12[4,:,:]),clev,cmap=cmap)
# 
# latplt = np.linspace(np.min(f12lat[5,:,:]), np.max(f12lat[5,:,:]), num=npx)
# lonplt = np.linspace(np.min(f12lon[5,:,:]), np.max(f12lon[5,:,:]), num=npx)
# lonpltm, latpltm = np.meshgrid(lonplt, latplt)
# x, y = map2(lonpltm,latpltm)
# cs2 = map2.contourf(x,y,np.rot90(f12[5,:,:]),clev,cmap=cmap)
# 
# latplt = np.linspace(np.min(f12lat[6,:,:]), np.max(f12lat[6,:,:]), num=npx)
# lonplt = np.linspace(np.min(f12lon[6,:,:]), np.max(f12lon[6,:,:]), num=npx)
# lonpltm, latpltm = np.meshgrid(lonplt, latplt)
# x, y = map3(lonpltm,latpltm)
# cs3 = map3.contourf(x,y,np.rot90(f12[6,:,:]),clev,cmap=cmap)
# 
# latplt = np.linspace(np.min(f12lat[7,:,:]), np.max(f12lat[7,:,:]), num=npx)
# lonplt = np.linspace(np.min(f12lon[7,:,:]), np.max(f12lon[7,:,:]), num=npx)
# lonpltm, latpltm = np.meshgrid(lonplt, latplt)
# x, y = map4(lonpltm,latpltm)
# cs4 = map4.contourf(x,y,np.rot90(f12[7,:,:]),clev,cmap=cmap)
# 
# fig.patch.set_facecolor('grey')
# 
# facemult = 1.115
# facedown = .016
# 
# plt.tight_layout()
# 
# pos1 = axN.get_position() 
# pos2 = [pos1.x0, pos1.y0,  pos1.width, pos1.height] 
# axN.set_position(pos2) 
# 
# pos1 = ax1.get_position() 
# pos2 = [pos1.x0+0.058, pos1.y0-facedown,  facemult*pos1.width, facemult*pos1.height] 
# ax1.set_position(pos2)
# 
# pos1 = ax2.get_position() 
# pos2 = [pos1.x0-0.014, pos1.y0-facedown,  facemult*pos1.width, facemult*pos1.height] 
# ax2.set_position(pos2)
# 
# pos1 = ax3.get_position() 
# pos2 = [pos1.x0-0.085, pos1.y0-facedown,  facemult*pos1.width, facemult*pos1.height] 
# ax3.set_position(pos2)
# 
# pos1 = ax4.get_position() 
# pos2 = [pos1.x0-0.155, pos1.y0-facedown,  facemult*pos1.width, facemult*pos1.height] 
# ax4.set_position(pos2)
# 
# pos1 = axS.get_position() 
# pos2 = [pos1.x0, pos1.y0,  pos1.width, pos1.height] 
# axS.set_position(pos2)
# 
# plt.savefig('CubedSpherePlotCoasts_Field-'+readvar+'_Level-'+str(plot_level)+'.'+filetype, bbox_inches='tight',facecolor=fig.get_facecolor())

#Nan out repeat locations
for tile in range(4):
    for j in range(npy):
        for i in range(npx):
            if i<j:
                f12[tile,j,i] = np.nan
            if i<-j+npx:
                f12[tile,j,i] = np.nan
            if i>j:
                f12[tile+8,j,i] = np.nan
            if i>-j+npx:
                f12[tile+8,j,i] = np.nan

#Cat to single array with nans
fp = np.zeros(12*npy*npx).reshape(4*npy,3*npx)
fp[:,0*npx+0:npx] = np.concatenate([f12[0,:,:], f12[1,:,:], f12[2 ,:,:], f12[3 ,:,:]])
fp[:,1*npx:2*npx] = np.concatenate([f12[4,:,:], f12[5,:,:], f12[6 ,:,:], f12[7 ,:,:]])
fp[:,2*npx:3*npx] = np.concatenate([f12[8,:,:], f12[9,:,:], f12[10,:,:], f12[11,:,:]])
fp = np.flipud(np.transpose(fp))


#Plot
fig = plt.figure(figsize=(14,8))
cp = plt.contourf(fp,clev,cmap=cmap)
cbar = plt.colorbar(cp,ticks=ctic)
if Dim2dor3d == '3D':
   plt.title('Cubed sphere plot of ' + readvar + itit + 'at level: ' + str(plot_level))
else:
   plt.title('Cubed sphere plot of ' + readvar)
plt.axis('off')
plt.axis('equal')

fig.patch.set_facecolor('grey')

plt.savefig('CubedSpherePlot_Field-'+readvar+'_Level-'+str(plot_level)+'.'+filetype, bbox_inches='tight',facecolor=fig.get_facecolor())
