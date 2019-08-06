#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 05:18:26 2019

  Plots using cartopy to show region of interest, grid scales, nesting, etc.

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

# bounds
ausbounds = [110,156,-44,-8] # EWSN
vicbounds = [135,154,-42,-33]
# projection type
proj=cartopy.crs.PlateCarree()

## Simple figure in first subplot:
##

# create figure, get current axes and set projection
fig = plt.figure(figsize=(12,12))

ax = fig.add_subplot(2, 2, 1, projection=proj)
#ax = plt.figure().gca(projection=cartopy.crs.PlateCarree())

# Add land, ocean, coasts, etc.
ax.add_feature(cpf.LAND)
ax.add_feature(cpf.OCEAN)
ax.add_feature(cpf.COASTLINE)
ax.add_feature(cpf.BORDERS, linestyle=':')
ax.add_feature(cpf.LAKES,   alpha=0.5)
ax.add_feature(cpf.RIVERS)

# Show just Australia [E,W, S,N]
ax.set_extent(ausbounds)

# title
plt.title("Plain")

## second subplot shows natural earth and text on top
##
ax = fig.add_subplot(2, 2, 2, projection=proj)

ax.set_extent([80, 170, -45, 30], crs=proj)

# Put a background image on for nice sea rendering.
ax.stock_img()

# Create a feature for States/Admin 1 regions at 1:50m from Natural Earth
states_provinces = cpf.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none')

SOURCE = 'Natural Earth'
LICENSE = 'public domain'

# add land, coastlines, states
ax.add_feature(cpf.LAND)
ax.add_feature(cpf.COASTLINE)
ax.add_feature(states_provinces, edgecolor='gray')

# Add a text annotation for the license information to the
# the bottom right corner.
text = matplotlib.offsetbox.AnchoredText(r'$\mathcircled{{c}}$ {}; license: {}'.format(SOURCE, LICENSE),
                    loc=4, prop={'size': 10}, frameon=True)
ax.add_artist(text)

# title
plt.title("Fancy")

## Shhow AUS with nested lat lon grid 
##


ax = fig.add_subplot(2, 1, 2, projection=proj)

# showing victoria
ax.set_extent(vicbounds, crs=proj)
plt.title("Victoria")

# Put a background image on for nice sea rendering.???
#ax.stock_img()


zoomedland = cartopy.feature.NaturalEarthFeature(category='physical',
                                                 name='land',
                                                 scale='10m',
                                                 facecolor="yellow")
zoomedcoast = cartopy.feature.NaturalEarthFeature(category='physical',
                                                  name='coastline',
                                                  scale='10m')
# levels: ist of integers 1-4 corresponding to the desired GSHHS feature levels to draw
zoomedcoast2 = cartopy.feature.GSHHSFeature(scale='auto',levels=[1,2,3,])

# add land, coastlines, states
#ax.add_feature(cpf.LAND)
ax.add_feature(zoomedland)
ax.add_feature(zoomedcoast2)
#ax.add_feature(cpf.COASTLINE)
# states?
#ax.add_feature(states_provinces, edgecolor='gray')


# 

# add grid at coarse resolution


# Create an inset GeoAxes showing the location of victoria. X,Y,width,heigh
sub_ax = plt.axes([0.151, 0.12, 0.15, 0.15], projection=proj)
sub_ax.set_extent(ausbounds, crs=proj)

# Make a nice border around the inset axes.
#effect = matplotlib.patheffects.Stroke(linewidth=4, foreground='wheat', alpha=0.5)
#sub_ax.outline_patch.set_path_effects([effect])

# Add the land, coastlines and the extent of victoria.
sub_ax.add_feature(cpf.LAND)
sub_ax.add_feature(cpf.COASTLINE)
# sgeom boxes based on ESWN?
extent_box = sgeom.box(vicbounds[0], vicbounds[2], vicbounds[1], vicbounds[3])
sub_ax.add_geometries([extent_box], proj, color='none',
                          edgecolor='blue', linewidth=2)

## finally show or save figure
##

plt.show()
