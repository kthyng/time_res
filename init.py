'''
Functions to initialize various numerical experiments.

Make a new init_* for your application.

loc     Path to directory of grid and output files
nsteps  Number of steps to do between model outputs (iter in tracmass)
ndays   number of days to track the particles from start date
ff      ff=1 to go forward in time and ff=-1 for backward in time
date    Start date in datetime object
tseas   Time between outputs in seconds
ah      Horizontal diffusion in m^2/s. 
        See project values of 350, 100, 0, 2000. For -turb,-diffusion
av      Vertical diffusion in m^2/s.
do3d    for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
doturb  turbulence/diffusion flag. 
        doturb=0 means no turb/diffusion,
        doturb=1 means adding parameterized turbulence
        doturb=2 means adding diffusion on a circle
        doturb=3 means adding diffusion on an ellipse (anisodiffusion)
lon0    Drifter starting locations in x/zonal direction.
lat0    Drifter starting locations in y/meridional direction.
z0/zpar Then z0 should be an array of initial drifter depths. 
        The array should be the same size as lon0 and be negative
        for under water. Currently drifter depths need to be above 
        the seabed for every x,y particle location for the script to run.
        To do 3D but start at surface, use z0=zeros(ia.shape) and have
         either zpar='fromMSL'
        choose fromMSL to have z0 starting depths be for that depth below the base 
        time-independent sea level (or mean sea level).
        choose 'fromZeta' to have z0 starting depths be for that depth below the
        time-dependent sea surface. Haven't quite finished the 'fromZeta' case.
        Then: 
        set z0 to 's' for 2D along a terrain-following slice
         and zpar to be the index of s level you want to use (0 to km-1)
        set z0 to 'rho' for 2D along a density surface
         and zpar to be the density value you want to use
         Can do the same thing with salinity ('salt') or temperature ('temp')
         The model output doesn't currently have density though.
        set z0 to 'z' for 2D along a depth slice
         and zpar to be the constant (negative) depth value you want to use
        To simulate drifters at the surface, set z0 to 's' 
         and zpar = grid['km']-1 to put them in the upper s level
         z0='s' is currently not working correctly!!!
         In the meantime, do surface using the 3d set up option but with 2d flag set
xp      x-locations in x,y coordinates for drifters
yp      y-locations in x,y coordinates for drifters
zp      z-locations (depths from mean sea level) for drifters
t       time for drifter tracks
name    Name of simulation to be used for netcdf file containing final tracks

'''

import numpy as np
import os
import netCDF4 as netCDF
import pdb
import glob
from datetime import datetime, timedelta
from matplotlib.mlab import *
import tracpy

units = 'seconds since 1970-01-01'

def init(grid):
    '''
    Initialization for seeding drifters at all shelf model grid points to be run
    forward.

    Optional inputs for making tests easy to run:
        date    Input date for name in datetime format
                e.g., datetime(2009, 11, 20, 0). If date not input,
                name will be 'temp' 
        grid    If input, will not redo this step. 
                Default is to load in grid.
    '''

    ndays = 25
    ff = 1 # forward in time
    ah = 0 # no diffusivity
    av = 0 # no diffusivity
    z0 = 's' # choosing by vertical level
    zpar = 59 # want top vertical level of grid
    zparuv = 0 # model output only has one level
    do3d = 0 # just 2d
    doturb = 0 # just comparing base paths
    dostream = 0 # don't do transport calculations
    loc = ['ocean_his_0001.nc', 'ocean_his_2010-07-01_00.nc']
    date = datetime(2010, 7, 1, 0, 5)
    # 5 minute resolution at highest, actual
    tseas = 5*60. # time between model outputs, in seconds
    units = 'seconds since 1970-01-01'    

    # Mesh of lat/lon starting points (same as horizontal_diffusivity, the LaCasce 2003 SCULP 1 locations)
    # llcrnrlon = -93.8; urcrnrlon = -92.2; llcrnrlat = 28; urcrnrlat = 29.2; # LaCasce (orig)
    llcrnrlon = -95.8; urcrnrlon = -94.2; llcrnrlat = 21.1; urcrnrlat = 22.3; # (4th)
    # llcrnrlon = -87.8; urcrnrlon = -86.2; llcrnrlat = 25; urcrnrlat = 26.2; # Test area
    xcrnrs, ycrnrs = grid['basemap']([llcrnrlon, urcrnrlon], [llcrnrlat, urcrnrlat])
    X, Y = np.meshgrid(np.arange(xcrnrs[0], xcrnrs[1], 2000), 
                        np.arange(ycrnrs[0], ycrnrs[1], 2000))
    # X, Y = np.meshgrid(np.arange(xcrnrs[0], xcrnrs[1], 700), 
    #                     np.arange(ycrnrs[0], ycrnrs[1], 700))
    lon0, lat0 = grid['basemap'](X, Y, inverse=True)

    # Eliminate points that are outside domain or in masked areas
    lon0, lat0 = tracpy.tools.check_points(lon0, lat0, grid)

    return ndays, ff, tseas, ah, av, lon0, lat0, \
            z0, zpar, zparuv, do3d, doturb, grid, dostream
