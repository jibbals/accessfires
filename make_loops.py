#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  20 2019

  Script to make the plots shown in animators.ipynb
  
@author: jesse
"""

import matplotlib
matplotlib.use('Agg')# don't plot on screen, send straight to file
# this is for NCI display issue

# allow input arguments
import argparse
import iris # use iris for better fio

# math and date stuff
import numpy as np
from datetime import datetime,timedelta

# multiproc
from multiprocessing import Pool
from time import sleep

# local modules
from utilities import plotting, utils, fio
from finished_plots import winds_2panel, clouds_2panel

def waroona_cloud_loop(dtime):
    '''
    make an hours worth of clouds_2panel plots starting at argument dtime
    First get iris cubes from each of the data files we will read,
        subset the data before reading it to save ram and run faster
        also read fire front at matching times
    then send all the data to plotting method
    '''
    um_hour=datetime(dtime.year,dtime.month,dtime.day,dtime.hour)
    extentname='waroona'
    vectorskip=9
    
    
    # Constraints on dimensions (save ram and reading time)
    West,East,South,North = plotting._extents_['waroona']
    constr_z = iris.Constraint(model_level_number=lambda cell: cell < 120)
    constr_lons = iris.Constraint(longitude = lambda cell: West <= cell <= East )
    constr_lats = iris.Constraint(latitude = lambda cell: South <= cell <= North )
    
    # metres [z, lat, lon]
    zro, = iris.load('data/waroona/umnsaa_2016010515_mdl_ro1.nc', ['height_above_reference_ellipsoid' &
                                                                   constr_z & 
                                                                   constr_lats & 
                                                                   constr_lons])
    
    topog = fio.read_nc_iris('data/waroona/umnsaa_2016010515_slv.nc',
                             constraints = 'surface_altitude'  & 
                                         constr_lats & 
                                         constr_lons)[0]
    
    
    # Read the cubes
    slv,ro1,th1,th2 = fio.read_waroona_iris(dtime=um_hour, 
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
    
    # Get wind speed cube using hypotenuse of u,v (I think this is the first action that actually reads any file data)
    s = iris.analysis.maths.apply_ufunc(np.hypot,u,v) 
    
    sh,Ta = th1.extract(['specific_humidity','air_temperature'])
    
    rh = utils.relative_humidity_from_specific(sh.data,Ta.data)
    
    qc1,qc2 = th2.extract(['mass_fraction_of_cloud_ice_in_air','mass_fraction_of_cloud_liquid_water_in_air'])
    
    qc = qc1+qc2
    
    ## fire front
    # read 6 time steps:
    ff_dtimes = np.array([um_hour + timedelta(hours=x/60.) for x in range(10,61,10)])
    ffpath = um_hour.strftime('data/waroona_fire/firefront.CSIRO_24h.%Y%m%dT1500Z.nc')
    ff = fio.read_fire(ffpath,dtimes=ff_dtimes,cube=True)
    ff = ff.extract(constr_lats & constr_lons) # subset lats,lons
    
    # datetime of outputs
    tdim = p.coord('time')
    d0 = datetime.strptime(str(tdim.units),'hours since %Y-%m-%d %H:%M:00')
    timesteps = utils.date_from_gregorian(tdim.points, d0=d0)
    
    # also loop over different transects
    for i_transect in np.arange(1,6.5,1, dtype=int):
        for tstep in range(len(timesteps)):
            clouds_2panel(topog.data, s[tstep,0].data, u[tstep,0].data, v[tstep,0].data,
                          qc[tstep].data, rh[tstep].data, 
                          ff[tstep].data,
                          zro.data, lat, lon,
                          dtime=timesteps[tstep],
                          extentname=extentname,
                          vectorskip=vectorskip, 
                          transect=i_transect)
    

def waroona_wind_loop(dtime,ff=None):
    '''
    Loop over transects over waroona and make the wind outline figures
    '''
    # Waroona loop
    extentname='waroona'
    vectorskip=9

    # List of hours for which we have output.. 0515 - 0614
    #date_list = [datetime(2016,1,5,15) + timedelta(hours=x) for x in range(24)]
    # first hour is already done, subsequent hours have some missing fields!
    topography,latt,lont = fio.read_topog('data/waroona/topog.nc')
    
    
    # Read the files for this hour
    waroona_outputs = fio.read_waroona(dtime,th2=False) # don't need cloud stuff for this plot
    slv, ro1, th1, th2 = waroona_outputs
    
    # grab the outputs desired
    waroona={}
    waroona['firefront'] = ff
    waroona['wind_speed'] = ro1['wind_speed']
    waroona['x_wind_destaggered'] = ro1['x_wind_destaggered']
    waroona['y_wind_destaggered'] = ro1['y_wind_destaggered']
    waroona['topog'] = topography
    waroona['zth'] = th1['zth']
    waroona['latitude'] = slv['latitude']
    waroona['longitude'] = slv['longitude']
    waroona['upward_air_velocity'] = th1['upward_air_velocity']
    waroona['time'] = th1['time']
    
    
    # datetime of outputs
    timesteps = utils.date_from_gregorian(th1['time'])
    
    # also loop over different transects
    for i_transect in np.arange(1,6.5,1, dtype=int):
        for tstep in range(len(timesteps)):
            winds_2panel(waroona,tstep=tstep,
                         extentname=extentname,
                         vectorskip=vectorskip, 
                         transect=i_transect)
    
    # Save ram hopefully
    del waroona, slv, ro1, th1, th2, waroona_outputs

##################################################################
#################  MAIN ##########################################
##################################################################

if __name__=='__main__':
    
    ## Input arguments
    parser = argparse.ArgumentParser()
    
    # Action flags
    parser.add_argument("--winds", 
                action="store_true", # save as args.clouds 
                help="run wind outline panel plots")
    parser.add_argument("--clouds", 
                action="store_true", 
                help="run cloud outline panel plots")
    # arguments with inputs
    parser.add_argument("--quarterday",
                type=int,
                help="which six hours will be running? (1,2,3, or 4)")
    args = parser.parse_args()
    
    qd = args.quarterday - 1
    assert (qd < 4) and (qd > -1), "--quarterday argument must be between 1 and 4 inclusive"
    rundateinds=np.arange(qd*6,qd*6+6,1,dtype=int)
    
    # First 24 hours:
    day1_dtimes = np.array([datetime(2016,1,5,15) + timedelta(hours=x) for x in range(24)])
    rundates = day1_dtimes[rundateinds]
    
    if args.winds:
        ## READ FIRE FRONT:
        ff, fft, latff, lonff = fio.read_fire()
        ffdt = utils.date_from_gregorian(fft/3600.,d0=datetime(2016,1,5,15))
        # ff is every 10 minutes, day1times is hourly
        # For everything except the first minute, we can send FF data 
    
    
    for i,dtime in enumerate(rundates):
        print("INFO: Running ",dtime)
        
        
        if args.clouds:
            waroona_cloud_loop(dtime)
    
        if args.winds:
            
            if dtime>day1_dtimes[0]:
                # some offset between firefront and model output
                ffoffind = 5+qd*6*6+(i-1)*6
                ffhour = ff[ffoffind:ffoffind+6]
                #print("DEBUG:",ffoffind,ffoffind+6)
                #print("DEBUG:",dtime,ffdt[ffoffind])
                assert (dtime.hour == ffdt[ffoffind].hour) and ffdt[ffoffind].minute==0, "fire front datetime doesn't match that of model output"
                waroona_wind_loop(dtime,ffhour)
            else:
                waroona_wind_loop(dtime)
    
    #with Pool(processes=2) as pool:
        
        ## Send each datetime to the process pool
        #pool.map( waroona_wind_loop, day1_dtimes )
        #pool.map( waroona_cloud_loop, day1_dtimes )
        
