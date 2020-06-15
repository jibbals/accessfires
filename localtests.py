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
import iris.quickplot as qplt

# read tiff file
from osgeo import gdal, osr
import cartopy.crs as ccrs

from utilities import utils, plotting, fio, constants

from scipy import interpolate

fdtime = datetime(2016,1,7,3)
fdtime2 = datetime(2016,1,7,4)
mr='waroona_run3'

cubes=fio.read_model_run(mr,[fdtime,fdtime2], add_z=True,add_winds=True)
print (cubes)

# currently different conversion script has different datetime format:

print(cubes[0].coord('time'))
print(cubes[1].coord('time'))
#z, w, u, topog = cubes.extract(['z_th','upward_air_velocity','u','surface_altitude'])
#lats,lons = w.coord('latitude').points, w.coord('longitude').points
#print(z.shape)
   
