#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 15:14:27 2020
    Look at surfaces in 3d
        show surface level, and atmospheric stuff in a set of 3d axes
@author: jesse
"""


import matplotlib
#matplotlib.use('Agg',warn=False)

# plotting stuff
import matplotlib.pyplot as plt
import numpy as np
import warnings
from datetime import datetime

# 3d plotting
# isosurface plotting available in plotly
import plotly.graph_objects as go
from plotly.io import show
import plotly.io as pio

# local modules
from utilities import plotting, utils, fio

###
## GLOBALS
###
_sn_ = 'threedee'

###
## METHODS
###


## Get run data
mr='waroona_run2'
top=12000 # km height
hours=fio.model_outputs[mr]['filedates']
#hour=hours[0]
hour = datetime(2016,1,6,10)
extentname = mr.split('_')[0]
extent = plotting._extents_[extentname]
cubes = fio.read_model_run(mr,fdtime=hour,extent=extent, 
                           add_theta=True, add_topog=True)
ff, = fio.read_fire(mr,dtimes=[hour],extent=extent)

print(cubes)

th, qc = cubes.extract(['potential_temperature','qc'])
topog, = cubes.extract(['surface_altitude'])
levh  = qc.coord('level_height').points
topind = np.sum(levh<top)
qc0= qc[0,:topind,:,:].data.data
th0= th[0,:topind,:,:].data.data
levh = levh[:topind]

# these are level, lat, lon cubes
lat,lon = qc.coord('latitude').points, qc.coord('longitude').points
lev = np.arange(len(qc0[:,0,0]))


qc0_xyz = np.moveaxis(np.moveaxis(qc0,0,2),0,1)
th0_xyz = np.moveaxis(np.moveaxis(th0,0,2),0,1)
topog_xy = topog.data.data.T
ff0_xy = ff[0].data.data

print("cloud min, max, sum>min:",np.min(qc0_xyz), np.max(qc0_xyz), np.sum(qc0_xyz>0.01))
print("Theta min, max, sum>min:",np.min(th0_xyz), np.max(th0_xyz), np.sum(th0_xyz>0.01))

## 3d figure:

X,Y,Z = np.meshgrid(lon,lat,levh) 
# X Y Z are now [lat, lon, lev] for some reason
[X,Y,Z] = [ np.moveaxis(arr,0,1) for arr in [X,Y,Z]]


## Topography
topog_surf = go.Surface(z=topog_xy,
                   x=X[:,:,0],
                   y=Y[:,:,0],
                   colorscale='earth',
                   showscale=False, # remove colour bar
                   )

## Fire contour
fire_line = go.Contour(#z=topog_xy,
                       #x=X[:,:,0],
                       #y=Y[:,:,0],
                       z=ff0_xy,
                       #value=ff0_xy,
                       colorscale='Hot',
                       #contours=dict(
                       #   start=-1,
                       #   end=0,
                       #   size=2)
                       )

## topography with fire as colour
topog_fire = go.Surface(z=topog_xy,
                   x=X[:,:,0],
                   y=Y[:,:,0],
                   colorscale='Hot',
                   surfacecolor=ff0_xy,
                   showscale=False, # remove colour bar
                   )

## Clouds
cloud_iso = go.Isosurface(x=X.flatten(),
                          y=Y.flatten(),
                          z=Z.flatten(),
                          value=qc0_xyz.flatten(),
                          isomin=.1,
                          isomax=1,
                          surface_count=3,
                          opacity=.7,
                          caps=dict(x_show=False,y_show=False))

# topog
fig0 = go.Figure(data=topog_surf)
show(fig0,renderer='browser')

# fire
figf = go.Figure(data=topog_fire)
show(figf,renderer='browser')

# clouds
fig1= go.Figure(data=cloud_iso)
show(fig1,renderer='browser')

# combined?
fig2 = go.Figure(data=[topog_fire,cloud_iso])
show(fig2,renderer='browser')

## Potential temperature
#fig2 = go.Figure(data=go.Isosurface(x=X.flatten(),
#                                    y=Y.flatten(),
#                                    z=Z.flatten(),
#                                    value=th0_xyz.flatten(),
#                                    isomin=310,
#                                    isomax=320,
#                                    #surface_fill=0.4,
#                                    caps=dict(x_show=False,y_show=False)))
#show(fig2, renderer='browser')


#if __name__ == '__main__':
#    plotting.init_plots()
#    mr = 'waroona_run1'
#    
#    hours=fio.model_outputs[mr]['filedates']
#    testhours = [hours[0]]
#    
#    make_plots_emberstorm(mr, hours=testhours)