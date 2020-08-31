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


mr='sirivan_run1'
umhours=fio.run_info[mr]['filedates'][:2]
cubes=fio.read_model_run(mr,extent=[149,150.2,-32.1,-31],fdtime=umhours[0],
                         HSkip=5)
print(cubes)
print(len(cubes[0].coord('latitude').points))
print(cubes[0].coord('latitude').points)
print(cubes[0].coord('latitude').points[::5])
#ff=fio.read_fire(mr)
#print(ff)


