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
import pandas

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

dt = datetime(2017,2,12,1)
sirivan_extent = [149.55, 149.56, -32.1, -32.09]
ff, = fio.read_fire('sirivan_run1',dtimes=[dt],extent=sirivan_extent)
cubes = fio.read_model_run('sirivan_run1',fdtime=dt, extent=sirivan_extent)
constraints = fio._constraints_from_extent_(sirivan_extent)
ff01, = fio.read_nc_iris('data/sirivan_run1/fire/fire_speed.01.nc',constraints=constraints)
print(ff01)
print(ff) # old run dx,dy = .004
print(cubes) # old run dx,dy = .004, 
lats = cubes[0].coord('latitude')
lons = cubes[0].coord('longitude')
flats = ff.coord('latitude')
flons = ff.coord('longitude')
flats1 = ff01.coord('latitude')
flons1 = ff01.coord('longitude')
for lat in [lats, flats, flats1]:
    print(lat)

for lon in [lons, flons, flons1]:
    print(lon)

## flats start one gridspace to the south, end together with lats
## flons start one gridspace to the east, end together with lons


assert np.all(np.isclose(lats.points, flats.points)), "old lats and flats don't match"
assert np.all(np.isclose(lons.points, flons.points)), "old lons and flons don't match"


