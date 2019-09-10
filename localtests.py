#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""

import matplotlib as mpl
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
from netCDF4 import Dataset

from datetime import datetime, timedelta
import timeit
import warnings

import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt
import matplotlib.patches as mpatches

from utilities import utils,fio,plotting




# Try reading file using iris
import iris
import iris.quickplot as qplt


#### SETUP PLOTTING:
plotting.init_plots()

### DATETIMES FOR UM AND FIRE OUTPUT
um_dtimes = np.array([datetime(2016,1,5,15,10) + timedelta(hours=x) for x in range(24)])

dtime0 = um_dtimes[0]
dtimeE = um_dtimes[-1]
#um_hour=datetime(dtime0.year,dtime0.month,dtime0.day,dtime0.hour)
um_hour=datetime(dtimeE.year,dtimeE.month,dtimeE.day,dtimeE.hour)
ff_dtimes = np.array([um_hour + timedelta(hours=x/60.) for x in range(10,61,10)])


##### TEST CODE HERE

# Read pc file using iris

#dtime, constraints=None, extent=None, add_winds=False, add_theta=False):
dtime=datetime(2016,1,6,8)
# read some data
cubes=fio.read_waroona_pcfile(dtime,extent=plotting._extents_['waroona'], add_winds=True,add_dewpoint=True)
print(cubes)


slv,ro1,th1,th2 = fio.read_waroona(dtime,extent=plotting._extents_['waroona'], add_dewpoint=True, add_winds=True)
press,temps,tempd  = th1.extract(['air_pressure','air_temperature','dewpoint_temperature'])
pro,u,v = ro1.extract(['air_pressure','u','v'])

## first interpolate pressure and temperature to latlon
lons=press.coord('longitude')
lats=press.coord('latitude')
press.convert_units('hPa') # convert to hPa
temps.convert_units('K') # convert to sir Kelvin
tempd.convert_units('K') # convert to kelvin







