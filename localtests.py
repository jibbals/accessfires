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

extent = plotting._extents_['waroona']
model_version='waroona_old'

dpath = fio.model_outputs[model_version]['path']
topogpath = fio.model_outputs[model_version]['topog']

allcubes=None

dtimes = [datetime(2016,1,6,7), datetime(2016,1,6,8)]
subdtimes = [datetime(2016,1,6,7,30) + timedelta(minutes=30*d) for d in range(3)]

#run1 = fio.read_waroona_run1(dtimes[0],extent=extent)
#print(run1)

print(dtimes)
print(subdtimes)
#subdtimes=None


run1 = fio.read_model_run('waroona_run1',fdtime=dtimes, subdtimes=subdtimes,
                          extent=extent, add_topog=True,
                          add_z=True,add_winds=True,add_theta=True,
                          add_dewpoint=True)
print(run1)
grabbed_times = utils.dates_from_iris(run1[0])
print(grabbed_times)

old = fio.read_model_run('waroona_old',fdtime=dtimes, subdtimes=subdtimes,
                         extent=extent, add_topog=True,
                         add_z=True,add_winds=True,add_theta=True,
                         add_dewpoint=True)
    
print(old)

grabbed_times = utils.dates_from_iris(old[0])
print(grabbed_times)
