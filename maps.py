# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 14:56:00 2019
    Script dealing with map figures for publication
@author: jgreensl
"""

# plotting 
from matplotlib import colors, ticker, patches
import matplotlib.pyplot as plt
# maths
import numpy as np

# mapping library
import cartopy

# geometry to add onto cartopy
import cartopy.feature
import shapely.geometry
import cartopy.crs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.transforms import offset_copy

# LOCAL IMPORTS
from utilities import utils, fio, constants, plotting

# Plot defaults
plotting.init_plots()


#f, ax, gproj = plotting.map_google(plotting._extents_['waroona'], zoom=11)

# Australia:
aust = [110,150,-40,-10]
# Request map from google
request = cimgt.GoogleTiles()
gproj=request.crs
# Use projection to set up plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=gproj)

# Where are we looking
ax.set_extent(aust)

# default interpolation ruins location names
ax.add_image(request, 3, interpolation='spline36') 

## ADD nested outlines

## Show grid resolution maybe with annotations
