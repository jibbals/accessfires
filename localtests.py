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
from matplotlib import patches, colors
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

fdtime = datetime(2016,1,7,3)
fdtime2 = datetime(2016,1,7,4)
mr='waroona_run3'

#cubes=fio.read_model_run(mr,[fdtime,fdtime2], add_z=True,add_winds=True)
#print (cubes)
#z, w, u, topog = cubes.extract(['z_th','upward_air_velocity','u','surface_altitude'])
#lats,lons = w.coord('latitude').points, w.coord('longitude').points
#print(z.shape)


from scipy.stats import norm
import matplotlib

# Select a color map
cmap = matplotlib.cm.bwr

# Some Test data
npts = 100
x = np.linspace(-4, 4, npts)
y = norm.pdf(x)
z = np.sin(2 * x)
normalize = colors.Normalize(vmin=z.min(), vmax=z.max())

# The plot
fig = plt.figure()
ax = fig.add_axes([0.12, 0.12, 0.68, 0.78])
plt.plot(x, y, color="gray")
for i in range(npts - 1):
    plt.fill_between([x[i], x[i+1]], [y[i], y[i+1]], color=cmap(normalize(z[i])))

cbax = fig.add_axes([0.85, 0.12, 0.05, 0.78])
cb = matplotlib.colorbar.ColorbarBase(cbax, cmap=cmap, norm=normalize, orientation='vertical')
cb.set_label("Sin function", rotation=270, labelpad=15)
plt.show()
   
