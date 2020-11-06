#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""


### IMPORTS

import numpy as np
import xarray as xr
import pandas as pd
import os, shutil
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import iris # interpolation etc
import time # check performance
from datetime import datetime
from utilities import utils, fio, constants, plotting


### GLOBALS
_sn_='localtests'



### CODE
run1="waroona_run3"
dt=datetime(2016,1,5,15)
extent=constants.extents['waroonaz']
cubes = fio.read_model_run(run1, fdtime=[dt], extent=extent, 
                            add_topog=True)                
topog=cubes.extract("surface_altitude")[0]
lats=topog.coord('latitude').points
lons=topog.coord('longitude').points
plotting.map_topography(extent,topog.data,lats,lons,cbar=False)
fio.save_fig_to_path("test.png",plt,)
