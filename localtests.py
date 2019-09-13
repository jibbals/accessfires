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

# read some data old old run Uncoupled
extent=plotting._extents_['waroona']
West,East,South,North = extent
constr_lons = iris.Constraint(longitude = lambda cell: West <= cell <= East )
constr_lats = iris.Constraint(latitude = lambda cell: South <= cell <= North )
constraints = constr_lats & constr_lons

orog_path = 'data/waroona_oldold/stage3_sfc_orog.nc'
orogcubes = fio.read_nc_iris(orog_path)#,constraints=constraints)
print(orogcubes)
#qplt.pcolormesh(orogcubes[0][0,0], vmax=500)
perth = slice(400,475,None), slice(400,425,None)
qplt.contourf(orogcubes[0][0,0][perth],40,vmin=-30,vmax=400, cmap='terrain')

xwind_path = 'data/waroona_oldold/combined_alltimes_xwind_regridded_stage3.nc'
xwindcubes = fio.read_nc_iris(xwind_path)#,constraints=constraints)
print(xwindcubes)
xwind = xwindcubes[0]
qplt.contourf(xwind[0,0][perth])

ywind_path = 'data/waroona_oldold/combined_alltimes_ywind_regridded_stage3.nc'
ywindcubes = fio.read_nc_iris(ywind_path,constraints=constraints)
print(ywindcubes)
ywind = ywindcubes[0]
qplt.contourf(ywind[0,0][perth])

