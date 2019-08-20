# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 08:58:26 2019
    Run tests on whatever in this file
@author: jgreensl
"""

import numpy as np
from matplotlib import pyplot as plt

from utilities import fio, utils, plotting


def compare_z_heights():
    dat = fio.read_pcfile('data/umnsaa_pc2016010515.nc')
    
    nlev=80 # Save on ram just look at lower levels
    tstep=0
    
    lat  = dat['latitude'   ][:]
    lon  = dat['longitude'  ][:]
    lat1 = dat['latitude_0' ][:] # Staggered
    lon1 = dat['longitude_0'][:] # Staggered
    z   = dat['model_level_number'  ][:]
    z1  = dat['model_level_number_0'][:]
    Ta  = dat['air_temperature'][:,:] # at surface
    p   = dat['air_pressure'][tstep,:nlev,:,:] # in pascals
    pmsl = dat['air_pressure_at_sea_level'][0,:,:] # just take surface 
    u1  = dat['x_wind'][tstep,:nlev,:,:] # wind speeds are on their directional grid edges (z,lats,lon1)
    v1  = dat['y_wind'][tstep,:nlev,:,:] # [z,lat1,lons]
    q   = dat['specific_humidity_0'][tstep,:nlev,:,:]
    lh   = dat['level_height'][:nlev] # in metres (grid box middle, distance from ground)
    w   = dat['upward_air_velocity'][tstep,:nlev,:,:]
    qc  = dat['mass_fraction_of_cloud_liquid_water_in_air'][tstep,:nlev,:,:] + dat['mass_fraction_of_cloud_ice_in_air'][tstep,:nlev,:,:]
    # Also some stuff based on calculated data (in fio.py)
    theta = dat['theta'][tstep,:nlev,:,:]
    zth = dat['zth'][tstep,:nlev,:,:]
    u   = dat['x_wind_destaggered'][tstep,:nlev,:,:]
    v   = dat['y_wind_destaggered'][tstep,:nlev,:,:]
    s   = dat['wind_speed'][tstep,:nlev,:,:]
    # Save some ram:
    del dat
    
    # Dimensions
    nz,ny,nx = p.shape
    

    # READ TOPOG DATA FROM PA
    topog, latt, lont = fio.read_topog()
    # 400x400 (not matching the pc file)
    
    # make actual height data
    replh = np.repeat(lh[:,np.newaxis],400,axis=1)
    replh = np.repeat(replh[:,:,np.newaxis],400,axis=2)
    h  = topog[np.newaxis,:,:] + replh
    
    # Compare height (m) to zth (m), should be similar
    # need to match lat/lon for proper look
    plt.scatter(h[:,50,50],zth[:,50,50], )
    plt.plot([0,h[-1,50,50]],[0,h[-1,50,50]],'k--',alpha=0.6,linewidth=2,label='1:1')
    plt.xlabel("model level height + topog")
    plt.ylabel("zth")
    plt.legend()
    plt.title("h vs zth, not index matched (yet)")