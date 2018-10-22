from netCDF4 import Dataset
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re

saveFigs = True

date = '20181004'
time = '000000'

fields_geos = ['ua','va','t','delp','q','qi','ql','o3mr']
fields_gfs  = ['ua','va','T','Ps','sphum','ice_wat','liq_wat','o3mr']

kind_gfs = ['fv_core','fv_core','fv_core','fv_core','fv_tracer','fv_tracer','fv_tracer','fv_tracer']

nm = len(fields_geos)

nl_type = 'geos'
nl1 = '../INPUTS/GEOS_c180/H21/GEOS.traj.c180.eta.'+date+'_'+time[0:4]+'z.nc4'
nl2 = '../INPUTS/GEOS_c180/H15/GEOS.traj.c180.eta.'+date+'_'+time[0:4]+'z.nc4'
   
tlm_type = 'gfs'
tlm = date+'.'+time+'.linearforecast.final.KIND.res.tileXX.nc'


if (nl_type == 'gfs'):
    tmp = re.sub("XX", str(1), nl1)
    fh_nl1 = Dataset(tmp, mode='r')
elif (nl_type == 'geos'):
    fh_nl1 = Dataset(nl1, mode='r')

# Read grid
# ---------
print('Reading grid')

if (nl_type == 'gfs'):
    nl_xid = 'xaxis_1'
    nl_yid = 'yaxis_1'
    nl_zid = 'zaxis_1'
elif (nl_type == 'geos'):
    nl_xid = 'lon'
    nl_yid = 'lat'
    nl_zid = 'lev'
 
im = len(fh_nl1.dimensions[nl_xid])
jm = len(fh_nl1.dimensions[nl_yid])
lm = len(fh_nl1.dimensions[nl_zid])

if (nl_type == 'gfs'):
    jm = jm * 6

jm_gfs = jm/6

print('Reading variables')

x_nlm = np.zeros(im*jm*lm*nm).reshape(nm,lm,jm,im)
x_tlm = np.zeros(im*jm*lm*nm).reshape(nm,lm,jm,im)

if (nl_type == 'gfs'):

    fields = fields_gfs

else:

    fields = fields_geos
    fh_nl2 = Dataset(nl2, mode='r')

    for i in range(nm):
        if (fields[i] != 'Ps'):
            x_nlm[i,:,:,:] = fh_nl2.variables[fields[i]][:] - fh_nl1.variables[fields[i]][:]
            if (fields[i] == 'delp'):
                tmp = np.sum(x_nlm[i,:,:,:],axis=0)
                for l in range(lm):
                    x_nlm[i,l,:,:] = tmp
        elif (fields[i] == 'Ps'):
            for l in range(lm):
                x_nlm[i,l,:,:] = fh_nl2.variables[fields[i]][:] - fh_nl1.variables[fields[i]][:]


if (tlm_type == 'gfs'):

    fields = fields_gfs

    for i in range(nm):

        for tile in range(6):

            filetile1 = re.sub("XX", str(tile+1), tlm)
            filetile = re.sub("KIND", kind_gfs[i], filetile1)

            fh_tlm = Dataset(filetile, mode = 'r')

            js = tile*jm_gfs
            je = (tile+1)*jm_gfs

            if (fields[i] != 'Ps'):
                x_tlm[i,:,js:je,:] = fh_tlm.variables[fields[i]][:]
                if (fields[i] == 'delp'):
                    tmp = np.sum(x_nlm[i,:,js:je,:],axis=0)
                    for l in range(lm):
                        x_tlm[i,l,js:je,:] = tmp
            elif (fields[i] == 'Ps'):
                for l in range(lm):
                    x_tlm[i,l,js:je,:] = fh_tlm.variables[fields[i]][:]


# Compute and plot correlations
# -----------------------------
print('Computing Correlations')
f, ax = plt.subplots(1, 8, figsize=(18, 6))

corr_tlm = np.zeros(lm*nm).reshape(nm,lm)
for i in range(nm):
    for j in range(lm):
        c1 = np.corrcoef(np.reshape(x_tlm[i,j,:,:],(im*jm,)),np.reshape(x_nlm[i,j,:,:],(im*jm,)))

        corr_tlm[i,j] = c1[0,1]

    rangeup1 = 0
    rangeup2 = 0
    for j in range(lm-1,-1,-1):
        if np.isnan(corr_tlm[i,j]) == True:
            rangeup1 = j+1
            break

    ax[i].plot(corr_tlm[i,rangeup1:lm],range(rangeup1+1,lm+1),color='b')
    ax[i].set_xlim(0, 1)
    ax[i].set_ylim(1, 72)
    ax[i].invert_yaxis()
    ax[i].set_title(fields_geos[i])
    if (i==0):
       ax[i].set_ylabel('Model Level')

    print('Corelation:',fields_geos[i],np.mean(corr_tlm[i,rangeup1:lm]) )

if saveFigs == True:
    plt.savefig('Correlation.png')

