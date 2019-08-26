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
from finished_plots import winds_2panel
from cloud_outline import waroona_cloud_loop
    

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
        ff, fft, latff, lonff = fio.read_fire_old()
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
        
