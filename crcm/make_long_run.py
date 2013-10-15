import netCDF4 as netCDF
import glob
import numpy as np
from datetime import datetime, timedelta
from matplotlib.mlab import find
import pdb
import op

units = 'seconds since 1970-01-01'
loc_output = '/pong/raid/data/crcm/BP/2010/' # crcm output

Files = glob.glob(loc_output + '07*/cpl_ocn*') # just take files from July
Files.sort()
Files = Files[:-3]
Files.insert(0,glob.glob(loc_output + '0629/cpl_ocn*')[0]) # so I can start right at July 1

# #TEMP
# Files = Files[1:3]

# Preprocess file to make names consistent, make times as numbers referenced to 1970-01-01, and have file name consistent.
# Also need to change the u and v locations to be on c grid. And input zeta as zero.
# Also, the stuff in on the coupling atmospheric grid, which has more points than the ocean grid. So, I need to exclude those points:
# Coupling/WRF grid: W-E: [1-856] x S-N: [1-811]
# Ocean grid: W-E: [117-771] x S-N: [22-510]

# Try making one file for each original, so concatenate here and do linear combination of 
# fields in time during overlap period.
for i, File in enumerate(Files[:-1]):

    d1 = netCDF.Dataset(File)
    d2 = netCDF.Dataset(Files[i+1])
    Times1 = d1.variables['Times'][:]
    Times2 = d2.variables['Times'][:]
    uin1 = d1.variables['OCN_US'][:, 21:510, 116:771];
    uin2 = d2.variables['OCN_US'][:, 21:510, 116:771];
    d1.close(); d2.close();

    # u and v are on rho grid to start
    yu = uin1.shape[1]; xu = uin1.shape[2] - 1; yv = uin1.shape[1] - 1; xv = uin1.shape[2]; nt = uin1.shape[0];
    zl = 1; # recreate all 60 levels so that tracpy code is consistent between model output and grid
    uin1 = op.resize(uin1, 2); # onto staggered grid
    uin1 = uin1.reshape((nt, 1, yu, xu)).repeat(zl,axis=1)
    uin2 = op.resize(uin2, 2); # onto staggered grid
    uin2 = uin2.reshape((nt, 1, yu, xu)).repeat(zl,axis=1)

    # Make own time vector
    date1 = datetime(2010, int(Times1[0,5] + Times1[0,6]), int(Times1[0,8] + Times1[0,9]), 
                    int(Times1[0,11] + Times1[0,12]), int(Times1[0,14] + Times1[0,15]))
    date2 = datetime(2010, int(Times2[0,5] + Times2[0,6]), int(Times2[0,8] + Times2[0,9]), 
                    int(Times2[0,11] + Times2[0,12]), int(Times2[0,14] + Times2[0,15]))
    t1 = np.zeros(nt); t2 = np.zeros(nt);
    for j in xrange(nt):
        t1[j] = netCDF.date2num(date1 + j*timedelta(minutes=5), units)
        t2[j] = netCDF.date2num(date2 + j*timedelta(minutes=5), units)
    dates1 = netCDF.num2date(t1, units)
    dates2 = netCDF.num2date(t2, units)

    # Find time indices of overlap period, 1 day
    # Make the second file a complete file
    # Indices for overlap period in dates1: istart to end of array
    istart_dates1 = find(dates1>=date2)[0]
    # Indices for overlap period in dates2: beginning of array to iend_dates2
    iend_dates2 = find(dates2>=dates1[-1])[0]+1 # +1 for proper array indexing
    # interpolation constant for each time step
    r2 = np.linspace(0, 1, t1[istart_dates1:].size)
    r1 = 1. - r2

    # Linearly combine model output in overlap region
    uin = uin2.copy()
    uin[:iend_dates2] = (r1.reshape(r1.size,1,1,1)*uin1[istart_dates1:] + r2.reshape(r1.size,1,1,1)*uin2[:iend_dates2])
    del(uin1,uin2)
    # pdb.set_trace()
    t = t2.copy()

    # Make new file
    rootgrp = netCDF.Dataset('ocean_his_' + str(i+1).zfill(4) + '.nc','w',format='NETCDF3_64BIT')
    rootgrp.createDimension('xu',xu)
    rootgrp.createDimension('yu',yu)
    rootgrp.createDimension('xv',xv)
    rootgrp.createDimension('yv',yv)
    rootgrp.createDimension('zl',zl)
    rootgrp.createDimension('nt',nt)
    ocean_time = rootgrp.createVariable('ocean_time','f8',('nt')) # 64-bit floating point
    u = rootgrp.createVariable('u','f4',('nt','zl','yu','xu')) # 64-bit floating point
    v = rootgrp.createVariable('v','f4',('nt','zl','yv','xv')) # 64-bit floating point
    u[:] = uin; ocean_time[:] = t;
    del(uin)

    d1 = netCDF.Dataset(File)
    d2 = netCDF.Dataset(Files[i+1])
    vin1 = d1.variables['OCN_VS'][:, 21:510, 116:771];
    vin2 = d2.variables['OCN_VS'][:, 21:510, 116:771];
    d1.close(); d2.close();
    # u and v are on rho grid to start
    vin1 = op.resize(vin1, 1); # onto staggered grid
    vin1 = vin1.reshape((nt, 1, yv, xv)).repeat(zl,axis=1)
    vin2 = op.resize(vin2, 1); # onto staggered grid
    vin2 = vin2.reshape((nt, 1, yv, xv)).repeat(zl,axis=1)

    # Linearly combine model output in overlap region
    vin = vin2.copy()
    vin[:iend_dates2] = (r1.reshape(r1.size,1,1,1)*vin1[istart_dates1:] + r2.reshape(r1.size,1,1,1)*vin2[:iend_dates2])
    del(vin1,vin2)

    v[:] = vin; del(vin)
    rootgrp.close()
