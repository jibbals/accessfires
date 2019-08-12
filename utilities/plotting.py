#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 14:06:49 2019

  In this script is repeatable generic plotting stuff like adding a scale to a cartopy map

@author: jesse
"""

import numpy as np
import cartopy
import matplotlib
from matplotlib import patheffects

###
### GLOBALS
###


# Extents: EWSN
_extents_             = {}
_extents_['waroona']  = [115.5,116.3, -33.1,-32.5] # local
_extents_['waroonas'] = [112,120,-35,-30] # synoptic
_extents_['yarloop']  = _extents_['waroona']
_extents_['yarloops'] = _extents_['waroonas']
_latlons_             = {}
_latlons_['waroona']  = -32.8430, 115.8526 # latlon of waroona
_latlons_['yarloop']  = -32.9534, 115.9124
_latlons_['perth']    = -31.9505, 115.8605


def init_plots():
  matplotlib.rcParams['font.size'] = 18.0
  matplotlib.rcParams["text.usetex"]      = False     # I forget what this is for, maybe allows latex labels?
  matplotlib.rcParams["legend.numpoints"] = 1         # one point for marker legends
  matplotlib.rcParams["figure.figsize"]   = (9, 7)    # Default figure size
  matplotlib.rcParams["font.size"]        = 17        # font sizes:
  matplotlib.rcParams["axes.titlesize"]   = 24        # title font size
  matplotlib.rcParams["axes.labelsize"]   = 19        #
  matplotlib.rcParams["xtick.labelsize"]  = 15        #
  matplotlib.rcParams["ytick.labelsize"]  = 15        #
  matplotlib.rcParams['image.cmap'] = 'plasma'        # Colormap default
  matplotlib.rcParams['axes.formatter.useoffset'] = False # another one I've forgotten the purpose of



def crass_section():
    '''
    Draw generic cross section based on 2D data, x axis, y axis, contours, cmap, norm, cbarformat
    '''
    
    
    # Contour inputs: xaxis, yaxis, data, colour gradient 
    plt.contourf(slicex,slicez,slicedata,slicelevels,cmap=cmap,norm=norm)
    plt.colorbar(format=cbarform)
    # Add contour lines
    plt.contour(slicex,slicez,slicedata,slicecontours,colors='k')            
    # make sure land is obvious
    plt.fill_between(xaxis,slicetopog,interpolate=True,facecolor='black')
    plt.xticks(None,None)
    if ztop != None:
        plt.ylim(0,ztop)
    
    
    
def utm_from_lon(lon):
    """
    utm_from_lon - UTM zone for a longitude

    Not right for some polar regions (Norway, Svalbard, Antartica)

    :param float lon: longitude
    :return: UTM zone number
    :rtype: int
    """
    return np.floor( ( lon + 180 ) / 6) + 1

def scale_bar(ax, proj, length, location=(0.5, 0.05), linewidth=3,
              units='km', m_per_unit=1000):
    """
    http://stackoverflow.com/a/35705477/1072212
    ax is the axes to draw the scalebar on.
    proj is the projection the axes are in
    location is center of the scalebar in axis coordinates ie. 0.5 is the middle of the plot
    length is the length of the scalebar in km.
    linewidth is the thickness of the scalebar.
    units is the name of the unit
    m_per_unit is the number of meters in a unit
    """
    # find lat/lon center to find best UTM zone
    #x0, x1, y0, y1 = ax.get_extent(proj.as_geodetic())
    x0, x1, y0, y1 = ax.get_extent(proj)
    # Projection in metres
    utm = cartopy.crs.UTM(utm_from_lon((x0+x1)/2))
    # Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(utm)
    # Turn the specified scalebar location into coordinates in metres
    sbcx, sbcy = x0 + (x1 - x0) * location[0], y0 + (y1 - y0) * location[1]
    # Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbcx - length * m_per_unit/2, sbcx + length * m_per_unit/2]
    # buffer for scalebar
    buffer = [patheffects.withStroke(linewidth=5, foreground="w")]
    # Plot the scalebar with buffer
    ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k',
        linewidth=linewidth, path_effects=buffer)
    # buffer for text
    buffer = [patheffects.withStroke(linewidth=3, foreground="w")]
    # Plot the scalebar label
    t0 = ax.text(sbcx, sbcy, str(length) + ' ' + units, transform=utm,
        horizontalalignment='center', verticalalignment='bottom',
        path_effects=buffer, zorder=2)
    left = x0+(x1-x0)*0.05
    # Plot the N arrow
    #t1 = ax.text(left, sbcy, u'\u25B2\nN', transform=utm,
    #    horizontalalignment='center', verticalalignment='bottom',
    #    path_effects=buffer, zorder=2)
    # Plot the scalebar without buffer, in case covered by text buffer
    ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k',
        linewidth=linewidth, zorder=3)
