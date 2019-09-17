#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker, colors, patches
import numpy as np
from netCDF4 import Dataset

from datetime import datetime, timedelta
import timeit
import warnings

import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt

from utilities import utils,fio,plotting


# Try reading file using iris
import iris
import iris.quickplot as qplt


#### SETUP PLOTTING:
plotting.init_plots()

### DATETIMES FOR UM AND FIRE OUTPUT
um_dtimes = np.array([datetime(2016,1,5,15,10) + timedelta(hours=x) for x in range(24)])
ff_dtimes = np.array([datetime(2016,1,5,15,1) +  timedelta(hours=x/60.) for x in range(1440)])
LT_offset = timedelta(hours=8) # AWST = UTC+8

dtime0 = um_dtimes[0]
dtimeE = um_dtimes[-1]



##### TEST CODE HERE

# Read fire output
extent = plotting._extents_['waroonaz'] # Super zoomed
frontal_exists_times = np.array([datetime(2016,1,5,15) + timedelta(hours=x) for x in range(10,24)])
FFront, SHeat, FSpeed = fio.read_fire(frontal_exists_times, extent=extent,
                                      firefront=True, sensibleheat=True, firespeed=True)
lon,lat = FFront.coord('longitude').points, FFront.coord('latitude').points
#print(FFront)
#print(SHeat)

# First plot topography
clists = fio.read_waroona(dtime0, extent=extent)
topog, = clists[0].extract('topog')

# fire contour colour map
cmap = mpl.cm.get_cmap('autumn')
# we have 14 hours with fire
rgba = cmap(np.linspace(0,1,14))

fig = plt.figure(figsize=[8,9])

ax1 = plt.subplot(2,1,1)
plt.contourf(lon,lat, topog.data, 50, cmap='terrain', vmin=-100)
plotting.map_add_locations(['waroona'],text=['Waroona'],dx=.03,dy=-.006)
# plot contours and label time stamp
for ii,dt in enumerate(frontal_exists_times):
    dstamp = (dt+LT_offset).strftime("%H(LT)")
    FF = FFront[ii].data
    assert np.min(FF) < 0, 'no fire yet'    
    cs = plt.contour(lon,lat,FFront[ii].data.T, np.array([0]), 
                     colors=[rgba[ii]],linewidths=1)
    # Add dstamp as inline label...
    if dt.hour%6 == 0:
        plt.clabel(cs, inline=True, fmt=dstamp, fontsize=9, fontcolor='k')


ax2 = plt.subplot(4,1,3)
maxflux=np.max(SHeat.data,axis=0)
cs=plt.contourf(lon, lat, maxflux.T, 30, 
                cmap='inferno', 
                locator=ticker.LogLocator())
plt.colorbar(cs, orientation='horizontal', pad=0)
plt.title('sensible heat flux (?)',y=.74)
plt.xticks([],[])
plt.yticks([],[])

ax3 = plt.subplot(4,1,4)
maxspeed=np.max(FSpeed.data,axis=0)
cs=plt.contourf(lon, lat, maxspeed.T, 30, 
                cmap=plotting._cmaps_['windspeed'])
plt.colorbar(cs, orientation='horizontal', pad=0)
plt.title('Max Firespeed (m/s?)', y=.74)
plt.xticks([],[])
plt.yticks([],[])


fio.save_fig('figures/waroona/fire_spread/waroona_run1.png',plt)