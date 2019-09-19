#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker, colors, patches
import numpy as np
from netCDF4 import Dataset

from datetime import datetime, timedelta
import timeit
import warnings

import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt

from utilities import utils,fio,plotting


# Try reading file using iris
import iris
import iris.quickplot as qplt
from iris.experimental.equalise_cubes import equalise_attributes

#### SETUP PLOTTING:
plotting.init_plots()

### DATETIMES FOR UM AND FIRE OUTPUT
um_dtimes = np.array([datetime(2016,1,5,15,10) + timedelta(hours=x) for x in range(24)])
ff_dtimes = np.array([datetime(2016,1,5,15,1) +  timedelta(hours=x/60.) for x in range(1440)])
LT_offset = timedelta(hours=8) # AWST = UTC+8

dtime0 = um_dtimes[0]
dtimeE = um_dtimes[-1]

##### TEST CODE HERE

allcubes=None
dtimes = [datetime(2016,1,6,6), datetime(2016,1,6,7)]
for dtime in dtimes:
    extent = plotting._extents_['waroona']
    
    cubelist = fio.read_waroona_pcfile(dtime, extent=extent, 
                                       add_zth=False, add_topog=False)
    #print(cubelist)
    if allcubes is None:
        allcubes = cubelist
    else:
        allcubes.extend(cubelist)
    
    # All the old run cubes have a time dim, so I can just concatenate each file

## First need to unify time dimension:
iris.util.unify_time_units(allcubes)
## Also need to equalise the attributes list
# I think this just deletes attributes which are not the same between matching cubes..
equalise_attributes(allcubes) 

## Join along the time dimension
allcubes = allcubes.concatenate()

## now subset to our subdtimes

## Now let's add z_th, and topography

## Add optional extras here
