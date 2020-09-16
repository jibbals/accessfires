#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""


# IMPORTS

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches, colors, cm
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection
from datetime import datetime, timedelta
import iris # so we can add level number constraint
import iris.analysis
import iris.quickplot as qplt

# read tiff file
from osgeo import gdal, osr
import cartopy.crs as ccrs

from utilities import utils, plotting, fio, constants

from scipy import interpolate

import pandas

_sn_='localtests'


mr='waroona_run3'
day1=fio.run_info[mr]['filedates'][:24]
cubes1=fio.read_model_run(mr,extent=plotting._extents_['waroona'],
                         fdtime=day1,
                         #HSkip=5,
                         )
print(cubes1)

w1 = cubes1.extract('upward_air_velocity')[0]
wmax1 = w1.collapsed(['time','model_level_number','longitude', 'latitude'], iris.analysis.MAX)

print("MAXIMUM:",wmax1.data)

day2=fio.run_info[mr]['filedates'][24:]
cubes2=fio.read_model_run(mr,extent=plotting._extents_['waroona'],
                         fdtime=day2,
                         #HSkip=5,
                         )
print(cubes2)

w2 = cubes2.extract('upward_air_velocity')[0]
wmax2 = w2.collapsed(['time','model_level_number','longitude', 'latitude'], iris.analysis.MAX)

print("MAXIMUM:",wmax2.data)
