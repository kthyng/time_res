import tracpy
import tracpy.calcs
import numpy
import glob
import netCDF4 as netCDF
import numpy as np

# Loop through tracks files and run calculations
Files = glob.glob('tracks/*.nc')
dc = netCDF.Dataset('tracks/tseas_use300_nsteps1.nc') # control case, 5 min output
lonpc = dc.variables['lonp'][:]; latpc = dc.variables['latp'][:]; tpc = dc.variables['tp'][:];

loc = ['crcm/', '/pong/raid/data/crcm/BP/2010/0703/ocean_his_0703_2010-07-03_00.nc']
grid = tracpy.inout.readgrid(loc)
for File in Files:

    # Read in tracks
    d = netCDF.Dataset(File)
    lonp = d.variables['lonp'][:]; latp = d.variables['latp'][:]; tp = d.variables['tp'][:];

    # change tracks to grid coords
    #loc = ['crcm/', '/pong/raid/data/crcm/BP/2010/0703/ocean_his_0703_2010-07-03_00.nc']
    #grid = tracpy.inout.readgrid(loc)
    tracpy.inout.save_ll2grid(File, grid)
    dg = netCDF.Dataset(File[:-3] + 'gc.nc')
    xp = dg.variables['xg'][:]; yp = dg.variables['yg'][:]; tp = dg.variables['tp'][:];

    # Numerical D
    D2, nnans = tracpy.calcs.rel_dispersion_comp(lonpc, latpc, tpc, lonp, latp, tp, r=1, squared=False)
    np.savez('calcs/' + File.split('/')[1][:-3] + 'Dcomp.npz', D2=D2, t=tp, nnans=nnans)

    # Numerical D^2
    D2, nnans = tracpy.calcs.rel_dispersion_comp(lonpc, latpc, tpc, lonp, latp, tp, r=1, squared=True)
    np.savez('calcs/' + File.split('/')[1][:-3] + 'D2comp.npz', D2=D2, t=tp, nnans=nnans)

    # Physical D
    D2, nnans = tracpy.calcs.rel_dispersion(lonp, latp, r=1, squared=False)
    np.savez('calcs/' + File.split('/')[1][:-3] + 'D.npz', D2=D2, t=tp, nnans=nnans)

    # Physical D^2
    D2, nnans = tracpy.calcs.rel_dispersion(lonp, latp, r=1, squared=True)
    np.savez('calcs/' + File.split('/')[1][:-3] + 'D2.npz', D2=D2, t=tp, nnans=nnans)

    # Physical a, absolute dispersion
    D2, nnans = tracpy.calcs.abs_dispersion(lonp, latp, squared=False)
    np.savez('calcs/' + File.split('/')[1][:-3] + 'a.npz', D2=D2, t=tp, nnans=nnans)

    # Physical a^2, absolute dispersion
    D2, nnans = tracpy.calcs.abs_dispersion(lonp, latp, squared=True)
    np.savez('calcs/' + File.split('/')[1][:-3] + 'a2.npz', D2=D2, t=tp, nnans=nnans)

    # Physical path length
    D2, nnans = tracpy.calcs.path(lonp, latp, squared=False)
    np.savez('calcs/' + File.split('/')[1][:-3] + 's.npz', D2=D2, t=tp, nnans=nnans)

    # Physical path length squared
    D2, nnans = tracpy.calcs.path(lonp, latp, squared=True)
    np.savez('calcs/' + File.split('/')[1][:-3] + 's2.npz', D2=D2, t=tp, nnans=nnans)

    # Calculate velocity along tracks
    varp = tracpy.calcs.Var(xp, yp, tp, 'u', dg)
    np.savez('calcs/' + File.split('/')[1][:-3] + 'u.npz', var=varp, t=tp, nnans=nnans)
    varp = tracpy.calcs.Var(xp, yp, tp, 'v', dg)
    np.savez('calcs/' + File.split('/')[1][:-3] + 'v.npz', var=varp, t=tp, nnans=nnans)
