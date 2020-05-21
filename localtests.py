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

cubes = iris.load("/g/data/en0/jwg574/iris_convert/converted_data/RA2M_astart.nc")

print(cubes)
keys_of_interest_3d = ["x_wind", "y_wind", "air_density",
                       "upward_air_velocity", "specific_humidity",
                       ]
keys_of_interest_2d = ["surface_temperature", "surface_altitude"]

for key in keys_of_interest_2d:
    keycubes = cubes.extract(key)
    print(keycubes)
    cube=keycubes[0]
    qplt.pcolormesh(cube)
    fio.save_fig("test","NYE",key,plt)
    
for key in keys_of_interest_3d:
    keycubes = cubes.extract(key)
    print(keycubes)
    cube=keycubes[0]
    levels=[0,5,10,25,50, 100]
    nlevs = len(levels)
    plt.figure(figsize=[12,12])
    for ilev, lev in enumerate(levels):
        plt.subplot(3, 2, ilev+1)
        qplt.pcolormesh(cube[lev])
        plt.title("model level %d"%lev)
    plt.suptitle(key)
    plt.tight_layout()
    fio.save_fig("test","NYE",key,plt)
    