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
import iris

from datetime import datetime, timedelta

from utilities import fio, plotting, utils

_sn_ = 'localtests'

import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from owslib.wmts import WebMapTileService

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy


## read tiff satellite image
#
from osgeo import gdal, osr

gdal.UseExceptions()


## First: read the geotiff image with GDAL.
fname = './data/Waroona_Landsat_8_2015.tiff'

ds = gdal.Open(fname)
data = ds.ReadAsArray()
gt = ds.GetGeoTransform()
proj = ds.GetProjection()

inproj = osr.SpatialReference()
inproj.ImportFromWkt(proj)

print("INFO: inproj")
print(inproj)

## Second: convert the WKT projection information to a cartopy projection.
projcs = inproj.GetAuthorityCode('PROJCS')
projection = ccrs.epsg(projcs)
print("INFO: projection")
print(projection)

## Third: the figure
subplot_kw = dict(projection=projection)
fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=subplot_kw)

extent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
          gt[3] + ds.RasterYSize * gt[5], gt[3])

img = ax.imshow(data[:3, :, :].transpose((1, 2, 0)), extent=extent,
                origin='upper')

## Create map projection
#ax = plt.axes(projection=ccrs.PlateCarree())
#ax.coastlines()

## put Aus in middle
#extent = plotting._extents_['waroonas']
#ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=135))
#ax.set_extent(extent)
#ax.gridlines()
#ax.coastlines(resolution='50m')




