# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 08:58:26 2019
    Run tests on whatever in this file
@author: jgreensl
"""

import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime, timedelta
import iris
import iris.quickplot as qplt
import timeit

from utilities import fio, utils, plotting

def read_time_comparison():
    '''
        does subsetting prior to reading save time? 
        YES
    '''
    #TEST FIRE STUFF
    ff_path = 'data/waroona_fire/firefront.csiro.01.nc'
    ff_dtimes = np.array([datetime(2016,1,5,15,1) + timedelta(hours=x/60.) for x in range(1440)])
    tstart = timeit.default_timer()
    ff = fio.read_fire(ff_path,dtimes=[ff_dtimes[1]], cube=True)
    print(ff.data.shape) # here is where actual read takes place
    read1 = timeit.default_timer() - tstart
    print("took %6.5f seconds to read %d timesteps"%(read1, len(ff.coord('time').points)))
    
    tstart = timeit.default_timer()
    ff = fio.read_fire(ff_path,dtimes=ff_dtimes[[1,11,25,36,46]], cube=True)
    print(ff.data.shape)
    read6 = timeit.default_timer() - tstart
    print("took %6.5f seconds to read %d timesteps"%(read6, len(ff.coord('time').points)))
    
    tstart = timeit.default_timer()
    ff = fio.read_fire(ff_path,dtimes=None, cube=True)
    print(ff.data.shape)
    read10th = timeit.default_timer() - tstart
    print("took %6.5f seconds to read %d timesteps"%(read10th, len(ff.coord('time').points)))

def check_height():
    '''
    look at z heights and see if they contain topography or not
    '''
    
    plotting.init_plots()
    
    um_dtimes = np.array([datetime(2016,1,5,15,10) + timedelta(hours=x) for x in range(24)])
    
    dtime = um_dtimes[0]
    um_hour=datetime(dtime.year,dtime.month,dtime.day,dtime.hour)
    
    # Constraints on dimensions
    West,East,South,North = plotting._extents_['waroona']
    constr_z = iris.Constraint(model_level_number=lambda cell: cell < 60)
    constr_lons = iris.Constraint(longitude = lambda cell: West <= cell <= East )
    constr_lats = iris.Constraint(latitude = lambda cell: South <= cell <= North )
    
    # metres [z, lat, lon]
    zro, = iris.load('data/waroona/umnsaa_2016010515_mdl_ro1.nc', ['height_above_reference_ellipsoid' &
                                                                   constr_z & 
                                                                   constr_lats & 
                                                                   constr_lons])
    
    topog, = fio.read_nc_iris('data/waroona/umnsaa_2016010515_slv.nc',
                             constraints = 'surface_altitude'  & 
                                         constr_lats & 
                                         constr_lons)
    
    plt.subplot(211)
    qplt.contourf(topog)
    plt.subplot(212)
    qplt.contourf(zro[0])

def check_interp():
    '''
    iris interpolation of wind speed fields
    '''
    
    um_dtimes = np.array([datetime(2016,1,5,15,10) + timedelta(hours=x) for x in range(24)])
    
    dtime = um_dtimes[0]
    # Constraints on dimensions
    West,East,South,North = plotting._extents_['waroona']
    constr_z = iris.Constraint(model_level_number=lambda cell: cell < 60)
    constr_lons = iris.Constraint(longitude = lambda cell: West <= cell <= East )
    constr_lats = iris.Constraint(latitude = lambda cell: South <= cell <= North )
    
    # Read the cubes
    slv,ro1,th1,th2 = fio.read_waroona_iris(dtime=dtime, 
                                            constraints = [constr_z &
                                                           constr_lons &
                                                           constr_lats])
    
    # wind speeds need to be interpolated onto non-staggered latlons
    # u is on lon1 dim
    # v is on lat1 dim
    
    # pull out wind speed
    p, u1, v1 = ro1.extract(['air_pressure','x_wind','y_wind'])
    # DESTAGGER u and v using iris interpolate
    # u1: [t,z, lat, lon1]
    # v1: [t,z, lat1, lon]  # put these both onto [t,z,lat,lon]
    u = u1.interpolate([('longitude',p.coord('longitude').points)],
                       iris.analysis.Linear())
    v = v1.interpolate([('latitude',p.coord('latitude').points)],
                       iris.analysis.Linear())
    lon=u.coord('longitude').points
    lat=u.coord('latitude').points
    lon1=u1.coord('longitude').points
    lat1=v1.coord('latitude').points
    plt.subplot(2,2,1)
    plt.contourf(lon1,lat,u1[0,0].data)
    plt.title('surface u1')
    
    plt.subplot(2,2,2)
    plt.contourf(lon,lat1,v1[0,0].data)
    plt.title('surface v1')
    
    plt.subplot(2,2,3)
    plt.contourf(lon,lat,u[0,0].data)
    plt.title('surface u')
    
    plt.subplot(2,2,4)
    plt.contourf(lon,lat,v[0,0].data)
    plt.title('surface v')
    


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