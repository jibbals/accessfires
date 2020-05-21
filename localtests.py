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
from matplotlib import patches
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import iris # so we can add level number constraint
import iris.analysis

# read tiff file
from osgeo import gdal, osr
import cartopy.crs as ccrs

from utilities import utils, plotting, fio, constants

model_run='waroona_run3'
extentname=model_run.split('_')[0]
# zoomed extent for analysis
extentnamef = extentname + 'f'
HSkip=None
extent=plotting._extents_[extentnamef]
localtime = False
proj_latlon = ccrs.PlateCarree() # lat,lon projection
qc_thresh = constants.cloud_threshold

dtimes = fio.model_outputs[model_run]['filedates']
u10,v10 = fio.read_fire(
    model_run=model_run, 
    extent=extent, 
    dtimes=dtimes, 
    firefront=False, 
    wind=True)

print(u10)
print(v10)
