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

from datetime import datetime

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


rpath = 'data/waroona/umnsaa_2016010515_mdl_ro1.nc'

rcubes = fio.read_nc_iris(rpath)
print(rcubes[0].name())
print(rcubes[0].units)

cubess = fio.read_waroona_iris()


fpath='data/umnsaa_pa2016010515.nc'
#cube = load_cubes(fpath)
cubes = iris.load(fpath, ['surface_altitude','land_binary_mask','latitude'])
cubes = iris.load(fpath, None)
print(cubes)
topog=cubes[11]
print(topog)
print(topog.data.shape)
print(topog.coord('latitude'))
print("INFO: Reading(iris) ",fpath, " ... ")

qplt.contourf(topog)
