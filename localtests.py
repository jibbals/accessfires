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

# read tiff file
from osgeo import gdal, osr
import cartopy.crs as ccrs

from utilities import utils, plotting, fio



## show tiff from qgis
#
#file='sirivan_map.tiff'

#show_grid=False
#locnames=None
#fig=None
#subplot_row_col_n=None
#subplot_extent=None
#""""""
#gdal.UseExceptions()
#
#path_to_tiff = "data/QGIS/"+file
#
#if subplot_row_col_n is None:
#    subplot_row_col_n = [1,1,1]
#
## gdal dataset
#ds = gdal.Open(path_to_tiff)
## RGBa image read in as a numpy array
#img = plt.imread(path_to_tiff)
#
## projection defined in QGIS
#projection = ccrs.epsg("3857")
#
## geotransform for tiff coords (?)
## tells us the image bounding coordinates
#gt = ds.GetGeoTransform()
#imgextent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
#             gt[3] + ds.RasterYSize * gt[5], gt[3])
#
#
#fig = plt.figure(figsize=(13, 9))
#ax = fig.add_subplot(111,projection=projection)
#ax.imshow(img,
#          extent=imgextent, 
#          origin='upper')
#plt.gca().set_extent(extent)
