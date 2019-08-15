#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 10:17:14 2019

    Plots to be animated

@author: jesse greenslade
"""

import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np

from utilities import plotting, utils

def winds_2panel(topog, w, s, u, v, z, lat, lon, extentname, transect=1, hourstamp = ""):
    '''
    211 Plot showing contour map and wind speed, along with near sites and transect
    223 plot showing vertical motion along transect
    224 plot showing wind speed along transect
    INPUTS:
        topography, vert motion, wind speed, z,lat,lon, extent
        transect = int from 1 to 6 for transect choice
    
    Since the boundaries will be set by the extent, I just plot all the locations in _latlons_ 
    '''
    # set font sizes and dpi
    plotting.init_plots()
    
    # get plot extent, and transect
    extent = plotting._extents_[extentname]
    start,end = plotting._transects_["%s%d"%(extentname,transect)]
    
    plt.figure(figsize=[7,10])
    ax1 = plt.subplot(2,1,1)
    
    # top panel is topography
    plotting.map_topography(extent,topog,lat,lon)
    plt.title('Topography, wind speed '+hourstamp)
    
    # start to end x=[lon0,lon1], y=[lat0, lat1]
    plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
             linewidth=2, marker='o', markersize=5)
    
    # add nearby towns
    plotting.map_add_locations(['waroona','yarloop','sirivan','uarbry'], 
                               text=['Waroona', 'Yarloop','Sir Ivan','Uarbry'], 
                               textcolor='lightgrey')
    # add fire ignition
    plotting.map_add_locations(['fire_waroona','fire_sirivan'],
                               text = ['Fire ignition']*2, 
                               color='r', marker='*', textcolor='lightgrey')
    # add pyroCB
    
    
    # Add vectors for winds
    # just surface, and one every 10 points to reduce density
    skip = (slice(None,None,10),slice(None,None,10))
    #mlon,mlat = np.meshgrid(lon,lat)
    plt.quiver(lon[skip[1]],lat[skip[0]],u[0][skip],v[0][skip], scale=60)
    
    ax2 = plt.subplot(2,2,3)
    plotting.transect_w(w,z, lat, lon,start,end,topog=topog)
    
    ax3 = plt.subplot(2,2,4)
    plotting.transect_s(s,zth,lat,lon,start,end,topog=topog,latt=latt,lont=lont)
    plt.yticks([])
    