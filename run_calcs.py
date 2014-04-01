import tracpy
import tracpy.calcs
import numpy
import glob
import netCDF4 as netCDF
import numpy as np
import pdb

locsave = '/Volumes/Emmons/projects/time_res/calcs/test/' #'calcs/' #'calcs/4th/' #test/'# '/Volumes/Emmons/projects/time_res/'

# Loop through tracks files and run calculations
floc = '/Volumes/Emmons/projects/time_res/tracks/test/' #4th/' #test/'
Files = glob.glob(floc + '*.nc')
# Files = glob.glob('tracks/*.nc')
dc = netCDF.Dataset('/Volumes/Emmons/projects/time_res/tracks/test/tseas_use300_nsteps1.nc') # control case, 5 min output
# dc = netCDF.Dataset('/Volumes/Emmons/projects/time_res/tracks/tseas_use300_nsteps1.nc') # control case, 5 min output
# dc = netCDF.Dataset('/Volumes/Emmons/projects/time_res/tracks/test/tseas_use300_nsteps1.nc') # control case, 5 min output
# # dc = netCDF.Dataset('tracks/tseas_use300_nsteps1.nc') # control case, 5 min output
lonpc = dc.variables['lonp'][:]; latpc = dc.variables['latp'][:]; tpc = dc.variables['tp'][:];

loc = ['/Volumes/Emmons/projects/time_res/crcm/', '/Volumes/Emmons/projects/time_res/forgrid/ocean_his_0703_2010-07-03_00.nc']
# loc = ['crcm/', '/pong/raid/data/crcm/BP/2010/0703/ocean_his_0703_2010-07-03_00.nc']
llcrnrlon=-97.944923400878906; urcrnrlon=-79.164894104003906;
llcrnrlat=18.132225036621094; urcrnrlat=30.847625732421861
grid = tracpy.inout.readgrid(loc, llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
                                  llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat)  # grid file, Gulf model domain
# grid = tracpy.inout.readgrid(loc[1], vert_filename=loc[0], llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
#                                   llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat)  # grid file, Gulf model domain
nc = netCDF.MFDataset('/Volumes/Emmons/projects/time_res/crcm/ocean*.nc')
dt = 1000 # number of time indices to loop through at once
# Files = [floc + 'tseas_use28800_nsteps96.nc', floc + 'tseas_use36000_nsteps120.nc']
# Files = [floc + 'tseas_use57600_nsteps192.nc', floc + 'tseas_use86400_nsteps288.nc']
for File in Files:

    if 'gc' in File:# or '10800' in File:
        continue

    # if '2700' not in File:
    #     continue

    # if '10800' not in File:
    #     continue
    print File
    # Read in tracks
    d = netCDF.Dataset(File)
    lonp = d.variables['lonp'][:]; latp = d.variables['latp'][:]; tp = d.variables['tp'][:];

    # Project into x and y coords
    xp, yp = grid['basemap'](lonp,latp)

    # # change tracks to grid coords
    # # tracpy.inout.save_ll2grid(floc, File.split('/')[-1], grid)
    # dg = netCDF.Dataset(File[:-3] + 'gc.nc')
    # xp = dg.variables['xg'][:]; yp = dg.variables['yg'][:]; tp = dg.variables['tp'][:];

    # # Calculate velocity along tracks
    # i = 0
    # up = np.ones(xp.shape)*np.nan
    # vp = np.ones(xp.shape)*np.nan
    # while i<=tp.size: # loop through time to save memory
    #     print i
    #     up[:,i:i+dt] = tracpy.calcs.Var(xp[:,i:i+dt], yp[:,i:i+dt], tp[i:i+dt], 'u', nc)
    #     vp[:,i:i+dt] = tracpy.calcs.Var(xp[:,i:i+dt], yp[:,i:i+dt], tp[i:i+dt], 'v', nc)
    #     i = i + dt
    # # pdb.set_trace()
    # print 'saving...'
    # np.savez(locsave + File.split('/')[-1][:-3] + 'u.npz', var=up, t=tp)
    # np.savez(locsave + File.split('/')[-1][:-3] + 'v.npz', var=vp, t=tp)
    # dg.close()

    # # Numerical D
    # D2, nnans, imax, imin, D2max, D2min = tracpy.calcs.rel_dispersion_comp(lonpc, latpc, 
    #                                         tpc, lonp, latp, tp, r=1, squared=False)
    # np.savez(locsave + File.split('/')[-1][:-3] + 'Dcomp.npz', D2=D2, t=tp, 
    #                         nnans=nnans, imax=imax, imin=imin, D2max=D2max, D2min=D2min)

    # Numerical D^2
    D2, nnans = tracpy.calcs.rel_dispersion_comp(lonpc, latpc, 
                                            tpc, lonp, latp, tp, r=1, squared=True)
    np.savez(locsave + File.split('/')[-1][:-3] + 'D2comp.npz', D2=D2, t=tp, nnans=nnans)
    # D2, nnans, imax, imin, D2max, D2min = tracpy.calcs.rel_dispersion_comp(lonpc, latpc, 
    #                                         tpc, lonp, latp, tp, r=1, squared=True)
    # np.savez(locsave + File.split('/')[-1][:-3] + 'D2comp.npz', D2=D2, t=tp, 
    #                         nnans=nnans, imax=imax, imin=imin, D2max=D2max, D2min=D2min)

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
    # np.savez(locsave + 'calcs/test/' + File.split('/')[-1][:-3] + 's.npz', D2=D2, t=tp, nnans=nnans)

    # # Physical path length squared
    # D2, nnans = tracpy.calcs.path(lonp, latp, squared=True)
    # np.savez(locsave + 'calcs/' + File.split('/')[-1][:-3] + 's2.npz', D2=D2, t=tp, nnans=nnans)

    # # Numerical D/s
    # D2, nnans, imax, imin, D2max, D2min = tracpy.calcs.rel_dispersion_comp_ds(lonpc, latpc,
    #                                          tpc, lonp, latp, tp, r=1, squared=False)
    # np.savez(locsave + File.split('/')[-1][:-3] + 'Dcompds.npz', D2=D2, t=tp, 
    #                             nnans=nnans, imax=imax, imin=imin, D2max=D2max, D2min=D2min)

    # # Numerical D/a
    # D2, nnans, imax, imin, D2max, D2min = tracpy.calcs.rel_dispersion_comp_da(lonpc, latpc,
    #                                          tpc, lonp, latp, tp, r=1, squared=False)
    # np.savez(locsave + File.split('/')[-1][:-3] + 'Dcompda.npz', D2=D2, t=tp, 
    #                             nnans=nnans, imax=imax, imin=imin, D2max=D2max, D2min=D2min)

    # # For a single drifter or a few drifters
    # # Numerical D/s
    # ind = np.array([624,814])
    # D2, nnans, imax, imin, D2max, D2min = tracpy.calcs.rel_dispersion_comp_ds(lonpc[ind,:], latpc[ind,:],
    #                                          tpc, lonp[ind,:], latp[ind,:], tp, r=1, squared=False)
    # np.savez(locsave + 'calcs/' + File.split('/')[-1][:-3] + 'Dcompds_624_814.npz', D2=D2, t=tp, 
    #                             nnans=nnans, imax=imax, imin=imin, D2max=D2max, D2min=D2min)

    # # Numerical D/a
    # D2, nnans, imax, imin, D2max, D2min = tracpy.calcs.rel_dispersion_comp_da(lonpc[ind,:], latpc[ind,:],
    #                                          tpc, lonp[ind,:], latp[ind,:], tp, r=1, squared=False)
    # np.savez(locsave + 'calcs/' + File.split('/')[-1][:-3] + 'Dcompda_624_814.npz', D2=D2, t=tp, 
    #                             nnans=nnans, imax=imax, imin=imin, D2max=D2max, D2min=D2min)

    # # 1st moments
    # M1x, nnans = tracpy.calcs.moment1(xp)
    # # np.savez(locsave + File.split('/')[-1][:-3] + 'M1x.npz', M=M1x, t=tp, nnans=nnans)
    # M1y, nnans = tracpy.calcs.moment1(yp)
    # # np.savez(locsave + File.split('/')[-1][:-3] + 'M1y.npz', M=M1y, t=tp, nnans=nnans)

    # 2nd moments
    # M2x, nnans = tracpy.calcs.moment2(xp, M1x)
    # np.savez(locsave + File.split('/')[-1][:-3] + 'M2x.npz', M=M2x, t=tp, nnans=nnans)
    # M2y, nnans = tracpy.calcs.moment2(yp, M1y)
    # np.savez(locsave + File.split('/')[-1][:-3] + 'M2y.npz', M=M2y, t=tp, nnans=nnans)

    # # 3rd moments
    # M3x = tracpy.calcs.moment3(xp, M1x)
    # np.savez(locsave + File.split('/')[-1][:-3] + 'M3x.npz', M=M3x, t=tp)
    # M3y = tracpy.calcs.moment3(yp, M1y)
    # np.savez(locsave + File.split('/')[-1][:-3] + 'M3y.npz', M=M3y, t=tp)
 
    # # 4th moments
    # M4x, nnans = tracpy.calcs.moment4(xp, M1x)
    # np.savez(locsave + File.split('/')[-1][:-3] + 'M4x.npz', M=M4x, t=tp, nnans=nnans)
    # M4y, nnans = tracpy.calcs.moment4(yp, M1y)
    # np.savez(locsave + File.split('/')[-1][:-3] + 'M4y.npz', M=M4y, t=tp, nnans=nnans)

    # # Trajectory skill score
    # l = min((lonpc.shape[0], lonp.shape[0]))
    # ss = tracpy.calcs.traj_ss(lonpc[:l,:], latpc[:l,:], lonp[:l,:], latp[:l,:])
    # np.savez(locsave + File.split('/')[-1][:-3] + 's.npz', ss=ss, t=tp)
    # # ssmean = ss.mean(axis=0)
    # # np.savez(locsave + File.split('/')[-1][:-3] + 'ssmean.npz', ssmean=ssmean, t=tp)

    # # Change from s to ss and get means
    # n = 0.1; nname = 'n0p1'
    # d = np.load(locsave + File.split('/')[-1][:-3] + 's.npz')
    # s = d['ss']; t=d['t']
    # smean = s.mean(axis=0)
    # # ind = (s>n)
    # # ss = 1-s/n
    # # ss[ind] = 0.
    # # ssmean = ss.mean(axis=0)
    # # np.savez(locsave + File.split('/')[-1][:-3] + 'ss_' + nname + '.npz', ss=ss, t=t, n=n)
    # # np.savez(locsave + File.split('/')[-1][:-3] + 'ssmean_' + nname + '.npz', ssmean=ssmean, t=t, n=n)
    # np.savez(locsave + File.split('/')[-1][:-3] + 'smean_' + nname + '.npz', smean=smean, t=t)
