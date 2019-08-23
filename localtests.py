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

import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt
import matplotlib.patches as mpatches

from utilities import utils,fio,plotting

# Try reading file using iris
import iris
import iris.quickplot as qplt
from iris.fileformats.netcdf import load_cubes

plotting.init_plots()

um_dtimes = np.array([datetime(2016,1,5,15,10) + timedelta(hours=x) for x in range(24)])

cubeslist = fio.read_waroona_iris(um_dtimes[0])

for cubes in cubeslist:
    print('========')
    print(cubes)


#fpath='C:/Users/jgreensl/Desktop/data/waroona_fire/firefront.CSIRO_24h.20160105T1500Z.nc'





