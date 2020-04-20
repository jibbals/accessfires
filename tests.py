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

def check_model_time_vs_temperature():
    """
    Make sure model time is UTC by checking temperature (should roughly match expected insolation)
    """
    # Check model temperature time series vs time stamp...

    r1 = fio.read_model_run('waroona_run1',add_topog=False)
    # get temperature
    T_surf = r1.extract('air_temperature')[0][:,0,:,:]
    
    dts = utils.dates_from_iris(T_surf) # UTC
    # WAST = UTC+8
    lts = [dts[i] + timedelta(hours=8) for i in range(len(dts))]
    print(T_surf)
    
    plt.close()
    plt.close()
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    for fig, ax, X, title in zip([fig1,fig2],[ax1,ax2],[dts,lts],['UTC','Local time']):
        plt.sca(ax)
        plt.plot_date(X, np.mean(T_surf.data,axis=(1,2)))
        
        # rotate and align the tick labels so they look better
        fig.autofmt_xdate()
    
        # use a more precise date string for the x axis locations
        #ax.fmt_xdata = dates.DateFormatter('%Y-%m-%d')
        plt.title('mean surface temperature '+title)

def wind_direction_calculation():
    """
    
    """
    # left, top, right, bot
    x= [1, 0, -1, 0]
    y= [0, 1, 0, -1]
    # calculate radians (-pi to pi range)
    rads = np.arctan2(y,x)
    assert np.all(np.isclose(rads, [0,np.pi/2.0, np.pi, -1*np.pi/2.0])), "rads doesn't make sense"
    # convert to anticlockwise from directly east
    wind_dir_math = (rads * 180/np.pi)%360
    assert np.all(np.isclose(wind_dir_math, [0, 90, 180, 270])), "math dir doesn't make sense"
    # meteorologicaly wind dir: 0 is due north, + is clockwise
    wind_dir = (-1 * wind_dir_math + 90 )%360
    assert np.all(np.isclose(wind_dir, [90, 0, 270, 180])), "met dir doesn't make sense"
    wind_dir2 = (-1*rads*180/np.pi+90)%360
    assert np.all(np.isclose(wind_dir2, [90, 0, 270, 180])), "met dir doesn't make sense"

def latlon_comparison_from_extent_subsetting():
    """
        fire and waroona_old runs were returning differently spaced arrays after subsetting
        fixed in fio._constraints_from_extent_()
    """
    for mr, extent, dt in zip(['waroona_old','sirivan_run1'], [plotting._extents_['waroona'],plotting._extents_['sirivan']], [datetime(2016,1,6,1),datetime(2017,2,12,1)]):
    
        ff, = fio.read_fire(mr,dtimes=[dt],extent=extent)
        cubes = fio.read_model_run(mr,fdtime=dt, extent=extent)
        print(ff) # old run dx,dy = .004
        print(cubes) # old run dx,dy = .004, 
        lats = cubes[0].coord('latitude').points
        lons = cubes[0].coord('longitude').points
        flats = ff.coord('latitude').points
        flons = ff.coord('longitude').points
        for lat in [lats, flats]:#, flats1]:
            print(lat[:2],'...',lat[-2:])
    
        for lon in [lons, flons]:#, flons1]:
            print(lon[:2],'...',lon[-2:])
        
        assert np.all(np.isclose(lats, flats)), "old lats and flats don't match"
        assert np.all(np.isclose(lons, flons)), "old lons and flons don't match"



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
    
def vorticity_test():
    """
    dummy wind fields, give expected vorticity?
    expected results was calced by hand
    """
    u = np.array([[1,1],[1,2],[2,3]])
    v = np.array([[1,2],[1,1],[2,3]])
    lats = np.array([3,2,1]) # rows
    lons = np.array([1,2]) # cols
    expected_zeta = np.array([[1,2],[.5,1],[2,2]])
    zeta, OW, OWZ = utils.vorticity(u,v,lats,lons)
    #print(zeta)
    assert np.all(expected_zeta == zeta), "Vorticity calculation does not match expectation"