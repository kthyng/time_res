'''
Script to run drifters backward from Galveston Bay to examine the Bay's 
connectivity with the shelf region.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
import init
from datetime import datetime, timedelta
import glob
from matplotlib.mlab import find

def run():

    # Location of preprocessed CRCM model output, on PONG
    loc = ['crcm/', '/pong/raid/data/crcm/BP/2010/0703/ocean_his_0703_2010-07-03_00.nc']

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('tracks/test'):
        os.makedirs('tracks/test')
    if not os.path.exists('figures'):
        os.makedirs('figures')
    if not os.path.exists('figures/test'):
        os.makedirs('figures/test')

    llcrnrlon=-97.944923400878906; urcrnrlon=-79.164894104003906;
    llcrnrlat=18.132225036621094; urcrnrlat=30.847625732421861
    grid = tracpy.inout.readgrid(loc, llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
                                      llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat)  # grid file, Gulf model domain

    date = np.array([datetime(2010, 7, 1, 0, 5)])

    # Set changing parameters to run through
    # time between model outputs, in seconds (5 min, 30 min, 1 hr, 2, 3, 4 hr, 6 hr)
    tseas_use = np.array([5, 10, 20, 30, 45, 60, 
                            60*2, 60*3, 60*4, 60*6])*60.
    # Number of linear interpolation steps in time (up to equivalent to 5 min output
    # so that the fields are updated at the same max times)
    nsteps = np.array([1, 2, 4, 6, 9, 12, 
                        24, 36, 48, 72]) 

    # Sampling timing. Make so that always sampling at 5 minutes, in all simulations.
    N = np.array([1, 2, 4, 6, 9, 12, 
                    24, 36, 48, 72]) 

    # loop through options
    for i in xrange(len(nsteps)):

        name = 'test/' + 'tseas_use' + str(int(tseas_use[i])) + '_nsteps' + str(nsteps[i]) # File names to use

        print 'simulation running: ' + name

        # If the particle trajectories have not been run, run them
        if not os.path.exists('tracks/' + name + '.nc'):

            print 'running tracks'

            # Read in simulation initialization
            ndays, ff, tseas, ah, av, lon0, lat0, z0, zpar, zparuv, do3d, doturb, \
                    grid, dostream = init.init(grid)
            lonp, latp, zp, t, grid = tracpy.run.run(loc, nsteps[i], ndays, ff, date, tseas, ah, av, lon0, lat0,
                                                     z0, zpar, do3d, doturb, name, grid=grid, dostream=dostream,
                                                     N=N[i], zparuv=zparuv, tseas_use=tseas_use[i])
        else:
            print 'skipping tracks'

        # If basic figures don't exist, make them
        if not os.path.exists('figures/' + name + 'tracks.png'):

            print 'running plots'

            # Read in and plot tracks
            d = netCDF.Dataset('tracks/' + name + '.nc')
            lonp = d.variables['lonp'][:]
            latp = d.variables['latp'][:]
            tracpy.plotting.tracks(lonp, latp, name, grid=grid)
            # tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')
            d.close()
        else:
            print 'skipping plots'
   


if __name__ == "__main__":
    run()    