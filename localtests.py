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
from iris.experimental.equalise_cubes import equalise_attributes

#### SETUP PLOTTING:
plotting.init_plots()

### DATETIMES FOR UM AND FIRE OUTPUT
um_dtimes = np.array([datetime(2016,1,5,15,10) + timedelta(hours=x) for x in range(24)])
ff_dtimes = np.array([datetime(2016,1,5,15,1) +  timedelta(hours=x/60.) for x in range(1440)])
LT_offset = timedelta(hours=8) # AWST = UTC+8

dtime0 = um_dtimes[0]
dtimeE = um_dtimes[-1]

##### TEST CODE HERE

extent = plotting._extents_['waroona']
model_version='waroona_old'

dpath = fio.model_outputs[model_version]['path']
topogpath = fio.model_outputs[model_version]['topog']

allcubes=None

dtimes = [datetime(2016,1,6,7), datetime(2016,1,6,8)]
subdtimes = [datetime(2016,1,6,7,30) + timedelta(minutes=30*d) for d in range(3)]
oldoldsubdtimes = [datetime(2016, 1, 6, 15, 15), datetime(2016, 1, 7, 6, 15), datetime(2016, 1, 7, 13, 45)]
#run1 = fio.read_waroona_run1(dtimes[0],extent=extent)
#print(run1)

print(dtimes)
print(subdtimes)
#subdtimes=None

cubes=fio.read_model_run('waroona_run1',fdtime=dtimes[0],extent=extent,add_z=True)
#print(cubes)



#### Look at height vs pressure
z,p=cubes.extract(['z_th','air_pressure'])
p.convert_units('hPa')
#p=np.mean(p.data,axis=0) # avg over time
# 140,lats,lons
z0 = np.mean(z.data,axis=(1,2)) # avg over lats,lons
p0 = np.mean(p.data,axis=(0,2,3)) # avg over time, lats, lons

checkslice=slice(50,70)
print(p0[checkslice])
print(z0[checkslice])



## WA EPSG = 4462
# URL of NASA GIBS: Global Imagery Browse Services

import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from owslib.wmts import WebMapTileService

import cartopy.crs as ccrs

URL = 'http://gibs.earthdata.nasa.gov/wmts/epsg4326/best/wmts.cgi'
wmts = WebMapTileService(URL)

# Layer for high res satellite image: https://wiki.earthdata.nasa.gov/display/GIBS/GIBS+Available+Imagery+Products#expand-SurfaceReflectance16Products
layer = 'Landsat_WELD_CorrectedReflectance_TrueColor_Global_Annual'

date_str = '2010-01-06' # landsat not available after 2010

## Plot setup
# plot projections (for coord ref transforms)
plot_CRS = ccrs.Mercator()
geodetic_CRS = ccrs.Geodetic()
# where are we looking?
extent = plotting._extents_['waroona'] # synoptic waroona
lon0, lon1, lat0, lat1 = extent
# transform to map corners
x0, y0 = plot_CRS.transform_point(lon0, lat0, geodetic_CRS)
x1, y1 = plot_CRS.transform_point(lon1, lat1, geodetic_CRS)
# keep aspect ratio based on lat lon corners
ysize = 8
xsize = ysize * (x1 - x0) / (y1 - y0)
fig = plt.figure(figsize=(xsize, ysize), dpi=100)
# create plot axes
ax = fig.add_axes([0, 0, 1, 1], projection=plot_CRS)
ax.set_xlim((x0, x1))
ax.set_ylim((y0, y1))
# add map tile from web service
ax.add_wmts(wmts, layer, wmts_kwargs={'time': date_str})
# add label for waroona
lat_w,lon_w = plotting._latlons_['waroona']
path_effects=[PathEffects.Stroke(linewidth=5, foreground='k'), PathEffects.Normal()]
ax.plot(lon_w, lat_w, color='grey', linewidth=0, marker='o', markersize=None,
        transform=geodetic_CRS,
        path_effects=path_effects)
txt = ax.text(lon_w, lat_w, 'Waroona', fontsize=8, color='wheat',
              transform=geodetic_CRS)
txt.set_path_effects([PathEffects.withStroke(linewidth=3,
                                             foreground='black')])
# add fire ignition
lat_f,lon_f = plotting._latlons_['fire_waroona']
cs=ax.plot(lon_f, lat_f, color='red', linewidth=0, marker='*', markersize=None,
           transform=geodetic_CRS,
           path_effects=path_effects)
# add layer name
lat_bl = lat0 + 0.02*(lat1-lat0)
lon_bl = lon0 + 0.05*(lon1-lon0)
txt = ax.text(lon_bl, lat_bl, wmts[layer].title, fontsize=14, color='wheat',
              transform=geodetic_CRS)
txt.set_path_effects([PathEffects.withStroke(linewidth=5,
                                             foreground='black')])

plt.show()


def working_gibsplot():
    URL = 'http://gibs.earthdata.nasa.gov/wmts/epsg4326/best/wmts.cgi'
    wmts = WebMapTileService(URL)
    
    # Layers for MODIS true color and snow RGB
    layers = ['MODIS_Terra_SurfaceReflectance_Bands143',
              'MODIS_Terra_CorrectedReflectance_Bands367']
    
    date_str = '2016-02-05'
    
    # Plot setup
    plot_CRS = ccrs.Mercator()
    geodetic_CRS = ccrs.Geodetic()
    x0, y0 = plot_CRS.transform_point(4.6, 43.1, geodetic_CRS)
    x1, y1 = plot_CRS.transform_point(11.0, 47.4, geodetic_CRS)
    ysize = 8
    xsize = 2 * ysize * (x1 - x0) / (y1 - y0)
    fig = plt.figure(figsize=(xsize, ysize), dpi=100)
    
    for layer, offset in zip(layers, [0, 0.5]):
        ax = fig.add_axes([offset, 0, 0.5, 1], projection=plot_CRS)
        ax.set_xlim((x0, x1))
        ax.set_ylim((y0, y1))
        ax.add_wmts(wmts, layer, wmts_kwargs={'time': date_str})
        txt = ax.text(4.7, 43.2, wmts[layer].title, fontsize=18, color='wheat',
                      transform=geodetic_CRS)
        txt.set_path_effects([PathEffects.withStroke(linewidth=5,
                                                     foreground='black')])
    plt.show()