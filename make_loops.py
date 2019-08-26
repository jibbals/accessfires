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

# math and date stuff
import numpy as np
from datetime import datetime,timedelta

# multiproc
from multiprocessing import Pool
from time import sleep

# local modules
from wind_outline import waroona_wind_loop
from cloud_outline import waroona_cloud_loop
    



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
    
    
    for i,dtime in enumerate(rundates):
        print("INFO: Running ",dtime)    
        
        if args.clouds:
            waroona_cloud_loop(dtime)
    
        if args.winds:
            waroona_wind_loop(dtime)
    
    #with Pool(processes=2) as pool:
        
        ## Send each datetime to the process pool
        #pool.map( waroona_wind_loop, day1_dtimes )
        #pool.map( waroona_cloud_loop, day1_dtimes )
        
