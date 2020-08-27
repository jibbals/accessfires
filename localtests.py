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


mr='sirivan_run6'
umhours=fio.run_info[mr]['filedates'][:2]
cubes=fio.read_model_run(mr,extent=[150.0,150.2,-31.99,-31.9],fdtime=umhours)
print(cubes)

ff=fio.read_fire(mr)
print(ff)


