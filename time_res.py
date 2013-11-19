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

# function to compute great circle distance between point lat1 and lon1 
# and arrays of points given by lons, lats or both same length arrays
# Haversine formula
def get_dist(lon1,lons,lat1,lats): 
    lon1 = lon1*np.pi/180.
    lons = lons*np.pi/180.
    lat1 = lat1*np.pi/180.
    lats = lats*np.pi/180.

    earth_radius = 6373.
    distance = earth_radius*2.0*np.arcsin(np.sqrt(np.sin(0.50*(lat1-lats))**2 \
                                       + np.cos(lat1)*np.cos(lats) \
                                       * np.sin(0.50*(lon1-lons))**2))
    return distance

def calc_dispersion(name, grid=None, r=1., ind=None, squared=False):
    '''
    This version is for comparing drifters starting at the same location over multiple 
    different simulations
    '''

    if grid is None:
        loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
        grid = tracpy.inout.readgrid(loc)
    else:
        grid = grid

    # Read in tracks to compare with control case
    d = netCDF.Dataset(name)
    # d = netCDF.Dataset('tracks/2006-02-01T00C_doturb2_ah20.nc')
    lonp = d.variables['lonp'][:]
    latp = d.variables['latp'][:]
    t = d.variables['tp'][:]
    d.close()

    # Read in tracks for control case (base)
    dbase = netCDF.Dataset('tracks/tseas_use300_nsteps1.nc')
    lonpbase = dbase.variables['lonp'][:]
    latpbase = dbase.variables['latp'][:]
    tbase = dbase.variables['tp'][:] # has extra entries
    dbase.close()

    # Account for time difference between two simulations
    tstride = int((t[1]-t[0])/(tbase[1]-tbase[0]))

    # Use tstride to equate arrays
    lonpbase = lonpbase[:,::tstride] # arrays end up different sizes due to tstride
    latpbase = latpbase[:,::tstride]
    tbase = tbase[::tstride] # has extra entries
    # pdb.set_trace()
    # base can be longer than others due to weirdness with stride in simulation
    lonpbase = lonpbase[:,:lonp.shape[1]]
    latpbase = latpbase[:,:lonp.shape[1]]
    tbase = tbase[:lonp.shape[1]] # has extra entries

    # Use input indices
    if ind is not None:
        lonpbase = lonpbase[ind,:]
        latpbase = latpbase[ind,:]
        lonp = lonp[ind,:]
        latp = latp[ind,:]

    # We know that drifters from the two sets have a one to one correspondence
    D2 = np.ones(lonp.shape[1])*np.nan
    nnans = np.zeros(lonp.shape[1]) # to collect number of non-nans over all drifters for a time
    for i in xrange(lonp.shape[0]): # loop through drifters, time is in array, axis=1
        # if i ==411:
        #     pdb.set_trace()
        dist = get_dist(lonpbase[i,:], lonp[i,:], 
                    latpbase[i,:], latp[i,:])
        if squared:
            D2 = np.nansum(np.vstack([D2, dist**2]), axis=0)
        else:
            D2 = np.nansum(np.vstack([D2, dist]), axis=0)
        nnans = nnans + ~np.isnan(dist)
    D2 = D2.squeeze()/nnans #len(pairs) # average over all pairs

    if squared:
        np.savez(name[:-3] + 'D2squared.npz', D2=D2, t=t, nnans=nnans)
    else:
        np.savez(name[:-3] + 'D2notsquared.npz', D2=D2, t=t, nnans=nnans)

    return D2, t


def calc_dispersion_onesim(name, grid=None, squared=False, r=1.):
    '''
    This version is for calculating dispersion amongst drifters in one simulation
    '''

    if grid is None:
        loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
        grid = tracpy.inout.readgrid(loc)
    else:
        grid = grid

    # Read in tracks
    d = netCDF.Dataset(name)
    # d = netCDF.Dataset('tracks/2006-02-01T00C_doturb2_ah20.nc')
    lonp = d.variables['lonp'][:]
    latp = d.variables['latp'][:]
    t = d.variables['tp'][:]
    d.close()

    dist = np.zeros((lonp.shape[0],lonp.shape[0]))
    for idrifter in xrange(lonp.shape[0]):
        # dist contains all of the distances from other drifters for each drifter
        dist[idrifter,:] = get_dist(lonp[idrifter,0], lonp[:,0], latp[idrifter,0], latp[:,0])

    # let the index in axis 0 be the drifter id
    ID = np.arange(lonp.shape[0])

    pairs = []
    for idrifter in xrange(lonp.shape[0]):
        ind = find(dist[idrifter,:]<=r)
        for i in ind:
            if ID[idrifter] != ID[i]:
                pairs.append([min(ID[idrifter], ID[i]), 
                                max(ID[idrifter], ID[i])])

    pairs_set = set(map(tuple,pairs))
    pairs = map(list,pairs_set)# now pairs has only unique pairs of drifters
    # pairs.sort() #unnecessary but handy for checking work

    D2 = np.ones(lonp.shape[1])*np.nan
    nnans = np.zeros(lonp.shape[1]) # to collect number of non-nans over all drifters for a time
    for ipair in xrange(len(pairs)):
        dist = get_dist(lonp[pairs[ipair][0],:], lonp[pairs[ipair][1],:], 
                    latp[pairs[ipair][0],:], latp[pairs[ipair][1],:])
        D2 = np.nansum(np.vstack([D2, dist**2]), axis=0)
        nnans = nnans + ~np.isnan(dist)
    D2 = D2.squeeze()/nnans #len(pairs) # average over all pairs

    return D2, t

def run_calc_dispersion(grid, squared=False):

    Files = glob.glob('tracks/*.nc')
    Files.sort()
    for File in Files:
        D2, t = calc_dispersion(File, grid=grid, squared=squared)


def run():

    # Location of preprocessed CRCM model output, on PONG
    loc = ['crcm/', '/pong/raid/data/crcm/BP/2010/0703/ocean_his_0703_2010-07-03_00.nc']

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('figures'):
        os.makedirs('figures')

    llcrnrlon=-97.944923400878906; urcrnrlon=-79.164894104003906;
    llcrnrlat=18.132225036621094; urcrnrlat=30.847625732421861
    grid = tracpy.inout.readgrid(loc, llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
                                      llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat)  # grid file, Gulf model domain

    date = np.array([datetime(2010, 7, 1, 0, 5)])

    # Set changing parameters to run through
    # time between model outputs, in seconds (5 min, 30 min, 1 hr, 4 hr, 6 hr, 8 hr)
    tseas_use = np.array([5, 10, 10, 20, 20, 30, 30, 60, 60, 
                            60*4, 60*4, 60*4, 60*4, 60*4, 60*4, 60*4, 60*4, 60*4, 
                            60*6, 60*6, 60*8])*60.
    # Number of linear interpolation steps in time (up to equivalent to 5 min output)
    nsteps = np.array([1, 1, 2, 1, 4, 1, 6, 1, 12, 
                            1, 2, 4, 6, 8, 16, 24, 36, 12*4, 
                            1, 12*6, 1]) 

    # loop through options
    for i in xrange(len(nsteps)):

        name = 'tseas_use' + str(int(tseas_use[i])) + '_nsteps' + str(nsteps[i]) # File names to use

        # If the particle trajectories have not been run, run them
        if not os.path.exists('tracks/' + name + '.nc'):

            # Read in simulation initialization
            ndays, ff, tseas, ah, av, lon0, lat0, z0, zpar, zparuv, do3d, doturb, \
                    grid, dostream = init.init(grid)
            lonp, latp, zp, t, grid = tracpy.run.run(loc, nsteps[i], ndays, ff, date, tseas, ah, av, lon0, lat0,
                                                     z0, zpar, do3d, doturb, name, grid=grid, dostream=dostream,
                                                     zparuv=zparuv, tseas_use=tseas_use[i])

        # If basic figures don't exist, make them
        if not os.path.exists('figures/' + name + '*.png'):

            # Read in and plot tracks
            d = netCDF.Dataset('tracks/' + name + '.nc')
            lonp = d.variables['lonp'][:]
            latp = d.variables['latp'][:]
            tracpy.plotting.tracks(lonp, latp, name, grid=grid)
            # tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')
            d.close()
   


if __name__ == "__main__":
    run()    