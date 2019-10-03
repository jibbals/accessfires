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

mr = 'waroona_old'
dt = datetime(2016,1,6,5)
ff, = fio.read_fire(mr,dtimes=[dt])
cubes = fio.read_model_run(mr,fdtime=dt)

print(ff)
print(cubes)
lats = cubes[0].coord('latitude').points
lons = cubes[0].coord('longitude').points
print(lats[:2],lats[-2:],lons[:2],lons[-2:])