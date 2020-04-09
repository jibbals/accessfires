#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""

#import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker, colors, patches
import numpy as np
import iris
from scipy.stats import gaussian_kde, cumfreq

from datetime import datetime, timedelta

from utilities import fio, plotting, utils

_sn_ = 'localtests'

import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from owslib.wmts import WebMapTileService

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy



import plotly.graph_objects as go
from plotly.io import show
import plotly.io as pio
# Make browser the default plotly renderer
pio.renderers.default = "browser"
## Surface heat colour scale for plotly
shcmax=20000.
shcmin=0
shcolor=[
    [0, 'rgb(255,255,255)'],  # start at white
    [100/shcmax, 'rgb(80,80,80)'], # head to grey for minor heat
    [1000/shcmax, 'rgb(255,0,0)'], # head to red for next milestone
    [10000/shcmax, 'rgb(255,0,255)'], # end at purple
    [1.0, 'rgb(0,0,0)'], # approach black off the scale I guess
    ]

import threedee

overfire=[116.09, 116.19, -32.91,-32.84]
model_run='waroona_run1'
hour=18
extent=overfire
extent=[115.88, 116.0, -32.86,-32.83]
HSkip=None
send_to_browser=True
"""
    look at waroona emberstorm bits
        surface wind speeds
        heat flux
        pot temp up to ~ 600m
        
"""
top_height=1000
hours=fio.model_outputs[model_run]['filedates']
dtime=hours[hour]
# [lat,lon], [lat1,lon1] of extent transect through middle
transect = [[(extent[2]+extent[3])/2,extent[0]],[(extent[2]+extent[3])/2,extent[1]]]

cubes = fio.read_model_run(model_run, 
                           fdtime=dtime, 
                           extent=extent, 
                           add_theta=True, 
                           add_topog=True, 
                           add_winds=True,
                           add_z=True,
                           HSkip=HSkip)

th, qc, z_th = cubes.extract(['potential_temperature','qc','z_th'])
# datetimes in hour output
cubetimes = utils.dates_from_iris(th)

ff, sh = fio.read_fire(model_run,
                       dtimes=cubetimes,
                       extent=extent,
                       sensibleheat=True,
                       HSkip=HSkip)

# Get the rest of the desired data
topog, = cubes.extract(['surface_altitude'])
d_topog = topog.data.data.T # convert to lon,lat
u,v = cubes.extract(['u','v'])
levh  = qc.coord('level_height').points
topind = np.sum(levh<top_height)

# these are level, lat, lon cubes
lat,lon = th.coord('latitude').points, th.coord('longitude').points

# dimensional mesh
X,Y,Z = np.meshgrid(lon,lat,levh) 
## X Y Z are now [lat, lon, lev] for some reason
[X,Y,Z] = [ np.moveaxis(arr,0,1) for arr in [X,Y,Z]]
## Now they are lon, lat, lev
## Cut down to desired level
[X, Y, Z] = [ arr[:,:,:topind] for arr in [X,Y,Z]]

for cubetime in cubetimes:
    surface_list=[]
    
    # surface sensible heat flux [t,lon,lat]
    d_sh = sh[0].data.data
    d_th = threedee.cube_to_xyz(th,ztopind=topind)
    # topography surface
    topog_layer = go.Surface(
        z=d_topog,
        x=X[:,:,0],
        y=Y[:,:,0],
        colorscale=shcolor, # surface heat colormap 
        cmin=shcmin,
        cmax=shcmax,
        #reversescale=True,
        surfacecolor=d_sh, # colour by sensible heat
        #opacityscale=[[0.0, 0], [100.0, .8], [shcmax, 1]], 
        #hidesurface=True,
        #showscale=False, # remove colour bar,
    )
    surface_list.append(topog_layer)
    
    ## Pull out transect to paint against the wall
    sliceth  = utils.cross_section(d_th,lat,lon,transect[0],transect[1],npoints=len(lon))
    slicetopog = utils.cross_section(d_topog,lat,lon,transect[0],transect[1],npoints=len(lon))
    slicez = utils.cross_section(z_th,lat,lon,transect[0],transect[1],npoints=len(lon))
    
    ### Paint figure
    figname = None
    if not send_to_browser:
        #figname = cubetime.strftime('figures/threedee/test_%Y%m%d%H%M.png')
        figname = fio.standard_fig_name(model_run,_sn_,cubetime,subdir='downslope')
    
    layoutargs = dict(title=cubetime.strftime('%Y%m%d%H%M(UTC)'),
                      font=dict(size=18, color="#111111"),)
                                
    create_figure(surface_list, filename=figname, **layoutargs)


