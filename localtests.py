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

# Try reading new outputs
dtime=datetime(2016,1,5,15)
waroona_outputs = fio.read_waroona(dtime)
slv, ro1, th1, th2 = waroona_outputs

# Check outputs
for i in range(4):
    for k,v in waroona_outputs[i].items():
        print("==== %d ===="%i)
        print(k, np.shape(v))