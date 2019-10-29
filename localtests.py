#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker, colors, patches, dates
import numpy as np


from datetime import datetime, timedelta

from utilities import fio, plotting, utils

_sn_ = 'localtests'


from cartopy import crs as ccrs
import cartopy.io.img_tiles as cimgt

linescan_extent = [149.48, 150.04, -32.18, -31.85]
subplot_extent2 = [0,0.02,1,0.45]
fig = plt.figure()
stamen_terrain = cimgt.Stamen('terrain-background')
ax2 = fig.add_axes(subplot_extent2, frameon=False, projection=stamen_terrain.crs)
ax2.set_extent(linescan_extent)
ax2.add_image(stamen_terrain, 10)