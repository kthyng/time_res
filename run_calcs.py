import tracpy
import tracpy.calcs
import numpy
import glob
import netCDF4 as netCDF
import numpy as np
import pdb

locsave = ''# '/Volumes/Emmons/projects/time_res/'

# Loop through tracks files and run calculations
Files = glob.glob('/Volumes/Emmons/projects/time_res/tracks/*.nc')
# Files = glob.glob('tracks/*.nc')
dc = netCDF.Dataset('/Volumes/Emmons/projects/time_res/tracks/tseas_use300_nsteps1.nc') # control case, 5 min output
# dc = netCDF.Dataset('tracks/tseas_use300_nsteps1.nc') # control case, 5 min output
lonpc = dc.variables['lonp'][:]; latpc = dc.variables['latp'][:]; tpc = dc.variables['tp'][:];

loc = ['/Volumes/Emmons/projects/time_res/crcm/', '/Volumes/Emmons/projects/time_res/forgrid/ocean_his_0703_2010-07-03_00.nc']
# loc = ['crcm/', '/pong/raid/data/crcm/BP/2010/0703/ocean_his_0703_2010-07-03_00.nc']
grid = tracpy.inout.readgrid(loc)

# nc = netCDF.MFDataset('/Volumes/Emmons/projects/time_res/crcm/ocean*.nc')
# dt = 1000 # number of time indices to loop through at once
for File in Files:

    if 'gc' in File:# or '10800' in File:
        continue

    if '10800' not in File:
        continue

    # Read in tracks
    d = netCDF.Dataset(File)
    lonp = d.variables['lonp'][:]; latp = d.variables['latp'][:]; tp = d.variables['tp'][:];

    # # change tracks to grid coords
    # # tracpy.inout.save_ll2grid(File, grid)
    # dg = netCDF.Dataset(File[:-3] + 'gc.nc')
    # xp = dg.variables['xg'][:]; yp = dg.variables['yg'][:]; tp = dg.variables['tp'][:];

    # # Calculate velocity along tracks
    # i = 0
    # up = np.ones(xp.shape)*np.nan
    # vp = np.ones(xp.shape)*np.nan
    # while i<=tp.size: # loop through time to save memory
    #     up[:,i:i+dt] = tracpy.calcs.Var(xp[:,i:i+dt], yp[:,i:i+dt], tp[i:i+dt], 'u', nc)
    #     vp[:,i:i+dt] = tracpy.calcs.Var(xp[:,i:i+dt], yp[:,i:i+dt], tp[i:i+dt], 'v', nc)
    #     i = i + dt
    # # pdb.set_trace()
    # np.savez(locsave + 'calcs/' + File.split('/')[-1][:-3] + 'u.npz', var=up, t=tp)
    # np.savez(locsave + 'calcs/' + File.split('/')[-1][:-3] + 'v.npz', var=vp, t=tp)
    # dg.close()

    # Numerical D
    D2, nnans, imax, imin, D2max, D2min = tracpy.calcs.rel_dispersion_comp(lonpc, latpc, 
                                            tpc, lonp, latp, tp, r=1, squared=False)
    np.savez(locsave + 'calcs/' + File.split('/')[-1][:-3] + 'Dcomp.npz', D2=D2, t=tp, 
                            nnans=nnans, imax=imax, imin=imin, D2max=D2max, D2min=D2min)

    # Numerical D^2
    D2, nnans, imax, imin, D2max, D2min = tracpy.calcs.rel_dispersion_comp(lonpc, latpc, 
                                            tpc, lonp, latp, tp, r=1, squared=True)
    np.savez(locsave + 'calcs/' + File.split('/')[-1][:-3] + 'D2comp.npz', D2=D2, t=tp, 
                            nnans=nnans, imax=imax, imin=imin, D2max=D2max, D2min=D2min)

    # # Physical D
    # D2, nnans = tracpy.calcs.rel_dispersion(lonp, latp, r=1, squared=False)
    # np.savez(locsave + 'calcs/' + File.split('/')[-1][:-3] + 'D.npz', D2=D2, t=tp, nnans=nnans)

    # # Physical D^2
    # D2, nnans = tracpy.calcs.rel_dispersion(lonp, latp, r=1, squared=True)
    # np.savez(locsave + 'calcs/' + File.split('/')[-1][:-3] + 'D2.npz', D2=D2, t=tp, nnans=nnans)

    # # Physical a, absolute dispersion
    # D2, nnans = tracpy.calcs.abs_dispersion(lonp, latp, squared=False)
    # np.savez(locsave + 'calcs/' + File.split('/')[-1][:-3] + 'a.npz', D2=D2, t=tp, nnans=nnans)

    # # Physical a^2, absolute dispersion
    # D2, nnans = tracpy.calcs.abs_dispersion(lonp, latp, squared=True)
    # np.savez(locsave + 'calcs/' + File.split('/')[-1][:-3] + 'a2.npz', D2=D2, t=tp, nnans=nnans)

    # # Physical path length
    # D2, nnans = tracpy.calcs.path(lonp, latp, squared=False)
    # np.savez(locsave + 'calcs/' + File.split('/')[-1][:-3] + 's.npz', D2=D2, t=tp, nnans=nnans)

    # # Physical path length squared
    # D2, nnans = tracpy.calcs.path(lonp, latp, squared=True)
    # np.savez(locsave + 'calcs/' + File.split('/')[-1][:-3] + 's2.npz', D2=D2, t=tp, nnans=nnans)

    # Numerical D/s
    D2, nnans, imax, imin, D2max, D2min = tracpy.calcs.rel_dispersion_comp_ds(lonpc, latpc,
                                             tpc, lonp, latp, tp, r=1, squared=False)
    np.savez(locsave + 'calcs/' + File.split('/')[-1][:-3] + 'Dcompds.npz', D2=D2, t=tp, 
                                nnans=nnans, imax=imax, imin=imin, D2max=D2max, D2min=D2min)

    # Numerical D/a
    D2, nnans, imax, imin, D2max, D2min = tracpy.calcs.rel_dispersion_comp_da(lonpc, latpc,
                                             tpc, lonp, latp, tp, r=1, squared=False)
    np.savez(locsave + 'calcs/' + File.split('/')[-1][:-3] + 'Dcompda.npz', D2=D2, t=tp, 
                                nnans=nnans, imax=imax, imin=imin, D2max=D2max, D2min=D2min)
