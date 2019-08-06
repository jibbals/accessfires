#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 10:17:14 2019

Quiver animation script

@author: jesse greenslade
"""


###
## IMPORT LIBRARIES
###

# cartopy for mapping
import cartopy
import cartopy.feature as cpf
from cartopy.io.img_tiles import Stamen # creative commons map images
# plotting library
import matplotlib.pyplot as plt
import matplotlib
# geometry to add onto cartopy
import shapely.geometry as sgeom
# numpy for arrays, math stuff
import numpy as np


###
## plot configuration 
###

# plot size
# title size
# label size




###
## DO STUFF
###


# Dummy wind vectors
lats,lons = np.linspace(-55,-25,20), np.linspace(125,155,30)
u = np.sin(np.meshgrid(lons,lats)[0]*np.pi/180.)
v = np.cos(np.meshgrid(lons,lats)[0]*np.pi/180.)
alter = np.random.normal(loc=0,scale=0.1,size=(2,20,30))
u2,v2 = u+alter[0], v+alter[1]
# Create map with terrain etc
ausextent = [110,156,-44,-8] # EWSN
vicextent = [135,154,-42,-33]
# projection type
#proj=cartopy.crs.PlateCarree()

# nice picture of ground/sea
#tiler = cartopy.io.img_tiles.MapQuestOpenAerial()
#proj = tiler.crs
proj = cartopy.crs.PlateCarree()
zoomedland = cartopy.feature.NaturalEarthFeature(category='physical',
                                                 name='land',
                                                 scale='10m',
                                                 facecolor="none")
zoomedcoast = cartopy.feature.NaturalEarthFeature(category='physical',
                                                  name='coastline',
                                                  scale='10m')


## top subplot shows wind speed quiver without/with coupling
fig = plt.figure(figsize=(11,11))
ax = fig.add_subplot(2, 1, 1, projection=proj)

# showing victoria
ax.set_extent(vicextent, crs=proj)
ax.add_feature(zoomedcoast)
ax.add_feature(zoomedland)
plt.title("Wind Vectors")

# Add arrows to show the wind vectors
plt.quiver(lons, lats, u, v, zorder=2,color='m', label='uncoupled')
plt.quiver(lons,lats, u2,v2, zorder=2,color='orange', label='coupled')
plt.legend(loc=1)

# 

# add grid at coarse resolution


# Create an inset GeoAxes showing the location of victoria. X,Y,width,heigh
sub_ax = plt.axes([0.157, 0.53, 0.15, 0.15], projection=proj)
sub_ax.set_extent(ausextent, crs=proj)

# Make a nice border around the inset axes.
#effect = matplotlib.patheffects.Stroke(linewidth=4, foreground='wheat', alpha=0.5)
#sub_ax.outline_patch.set_path_effects([effect])

# Add the land, coastlines and the extent of victoria.
sub_ax.add_feature(cpf.LAND)
sub_ax.add_feature(cpf.COASTLINE)
# sgeom boxes based on ESWN?
extent_box = sgeom.box(vicextent[0], vicextent[2], vicextent[1], vicextent[3])
sub_ax.add_geometries([extent_box], proj, color='none',
                          edgecolor='blue', linewidth=2)

## finally show or save figure
##

plt.show()