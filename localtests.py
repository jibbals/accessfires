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

model_run='waroona_run1'
extentname=model_run.split('_')[0]
# zoomed extent for analysis
extentnamef = extentname + 'f'
HSkip=None
extent=plotting._extents_[extentnamef]
localtime = False
proj_latlon = ccrs.PlateCarree() # lat,lon projection
qc_thresh = constants.cloud_threshold


# READ EVERYTHING, SUBSET, CALC DESIRED METRICS, PLOT METRICS
cubes = fio.read_model_run(model_run, extent=extent)
ctimes = utils.dates_from_iris(cubes[0])
    
## get temperature, RH, cloud
q,T,qc = cubes.extract(['specific_humidity','air_temperature','qc'])
    