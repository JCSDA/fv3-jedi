from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

#User input required for the follwing:
plot_diff = 1         #Plot path2/file - path1/file
path1  = '/discover/nobackup/drholdaw/Jedi/Experiments/Test12/'  #Path of first file
path2  = '/discover/nobackup/drholdaw/Jedi/Experiments/Test13/'  #Path of second file
file_tplt_befr = 'fv_core.res.tile'  #Filename before tile number
file_tplt_aftr = '.nc'               #Filename after tile number
xdimvar = 'xaxis_1'                  #What to read to get dimension
ydimvar = 'yaxis_2'                  #What to read to get dimension
zdimvar = 'zaxis_1'                  #What to read to get dimension
readvar = 'T'                        #Variable to plot
Dim2dor3d = '3D'                     #Is this 2D or 3D field?
plot_level = 40                      #If 3D plot this level

#Another example:
#path  = '/discover/nobackup/drholdaw/Jedi/Experiments/INPUT/'
#file_tplt_befr = 'C96_grid_spec.tile'
#xdimvar = 'grid_x'
#ydimvar = 'grid_y'
#zdimvar = ' '
#readvar = 'grid_lon'
#Dim2dor3d = '2D'

#Tile 1 handle for getting dimension
fh1 = Dataset(path1 + file_tplt_befr + str(1) + file_tplt_aftr, mode='r')

#Dimensions
npx = len(fh1.dimensions[xdimvar])
npy = len(fh1.dimensions[ydimvar])
npz = 0
if Dim2dor3d == '3D':
   npz = len(fh1.dimensions[zdimvar])
   fr = np.zeros(6*npz*npy*npx).reshape(6,npz,npy,npx)
else:
   fr = np.zeros(6*npy*npx).reshape(6,npy,npx)

#Read file
for tile in range(6):
    file_tile = file_tplt_befr + str(tile+1) + file_tplt_aftr
    pathfile1 = path1 + file_tile
    pathfile2 = path2 + file_tile
    print(pathfile1)
    if plot_diff == 1:
        print(' '+pathfile2)
    fh1 = Dataset(pathfile1, mode='r')
    fh2 = Dataset(pathfile2, mode='r')
    if plot_diff == 1:
        fr[tile,:,:] = fh2.variables[readvar][:] - fh1.variables[readvar][:]
    else:
        fr[tile,:,:] = fh1.variables[readvar][:]

#Get level to plot
f = np.zeros(6*npy*npx).reshape(6,npy,npx)
if Dim2dor3d == '3D':
   f[:,:,:] = fr[:,plot_level,:,:]
else:
   f = fr

#Initial arrays
f12 = np.zeros(12*npy*npx).reshape(12,npy,npx)
fp = np.zeros(12*npy*npx).reshape(4*npy,3*npx)

#Transpose, rotate etc
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

#Cat to single array
fp[:,0*npx+0:npx] = np.concatenate([f12[0,:,:], f12[1,:,:], f12[2 ,:,:], f12[3 ,:,:]])
fp[:,1*npx:2*npx] = np.concatenate([f12[4,:,:], f12[5,:,:], f12[6 ,:,:], f12[7 ,:,:]])
fp[:,2*npx:3*npx] = np.concatenate([f12[8,:,:], f12[9,:,:], f12[10,:,:], f12[11,:,:]])
fp = np.flipud(np.transpose(fp))

#Contour levels
maxf = np.nanmax(fp)
minf = np.nanmin(fp)

if minf < 0:
    minf = -maxf
incf = (maxf - minf)/51
clev = np.arange(minf,maxf+incf,incf)

#Plot
fig = plt.figure(figsize=(14,8))
cp = plt.contourf(fp,clev)
cbar = plt.colorbar(cp)
if Dim2dor3d == '3D':
   plt.title('Cubed sphere plot of ' + readvar + ' at level: ' + str(plot_level))
else:
   plt.title('Cubed sphere plot of ' + readvar)
plt.axis('off')
plt.axis('equal')

plt.savefig('CubedSpherePlot_Field-'+readvar+'_Level-'+str(plot_level)+'.png', bbox_inches='tight')
