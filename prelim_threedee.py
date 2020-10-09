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

from utilities import fio, plotting, utils, constants

_sn_ = 'localtests'

import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from owslib.wmts import WebMapTileService

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy


####
# three dim testing
import plotly.graph_objects as go
import plotly.io as pio
# Make browser the default plotly renderer
pio.renderers.default = "browser"
    
### Get run data
mr='waroona_run2'
top=8000 # cloud top height to be shown
top_th=1000 # theta top height 
hours=fio.model_outputs[mr]['filedates']
hour=hours[0]
hour = datetime(2016,1,6,10)
extentname = mr.split('_')[0]
extent = constants.extents[extentname]
northerly_view = dict(scene_camera_eye=dict(x=0,y=-1,z=1.5))
standard_view = dict(scene_camera_eye=dict(x=1.3,y=-1.3,z=.6))

cubes = fio.read_model_run(mr,fdtime=hour,extent=extent, 
                           add_theta=True, add_topog=True, add_winds=True)
ff, = fio.read_fire(mr,dtimes=[hour],extent=extent)

print(cubes)

th, qc = cubes.extract(['potential_temperature','qc'])
topog, = cubes.extract(['surface_altitude'])
u,v,w = cubes.extract(['u','v','upward_air_velocity'])
levh  = qc.coord('level_height').points
topind = np.sum(levh<top)
topind_th = np.sum(levh<top_th)

qc0= qc[0,:topind,:,:].data.data
th0= th[0,:topind,:,:].data.data
levh = levh[:topind]
u0= u[0,:topind,:,:].data.data
v0= v[0,:topind,:,:].data.data
w0= w[0,:topind,:,:].data.data

# these are level, lat, lon cubes
lat,lon = qc.coord('latitude').points, qc.coord('longitude').points
lev = np.arange(len(qc0[:,0,0]))


qc0_xyz = np.moveaxis(np.moveaxis(qc0,0,2),0,1)
th0_xyz = np.moveaxis(np.moveaxis(th0,0,2),0,1)
u0_xyz = np.moveaxis(np.moveaxis(u0,0,2),0,1)
v0_xyz = np.moveaxis(np.moveaxis(v0,0,2),0,1)
w0_xyz = np.moveaxis(np.moveaxis(w0,0,2),0,1)
topog_xy = topog.data.data.T
ff0_xy = ff[0].data.data

#
### 3d figure:
#
X,Y,Z = np.meshgrid(lon,lat,levh) 
#    
## X Y Z are now [lat, lon, lev] for some reason
[X,Y,Z] = [ np.moveaxis(arr,0,1) for arr in [X,Y,Z]]
## Now they are lon, lat, lev
#
## Add topography to Z
Zh=np.zeros(np.shape(Z))
for i in range(len(levh)):
    Zh[:,:,i] = Z[:,:,i]+topog_xy



### Topography
## surface given by topography 2D output\n",
topog_surf = go.Surface(
    z=topog_xy,
    x=X[:,:,0],
    y=Y[:,:,0],
    colorscale='earth_r',
    showscale=False, # remove colour bar,
    )

flat_topog = go.Surface(
    z=Z[:,:,0],
    x=X[:,:,0],
    y=Y[:,:,0],
    colorscale='earth',
    surfacecolor=topog_xy,
    showscale=False, # remove colour bar,
    )

## topography with fire as colour
topog_fire = go.Surface(
    z=topog_xy,
    x=X[:,:,0],
    y=Y[:,:,0],
    colorscale='Hot',
    surfacecolor=ff0_xy,
    showscale=False, # remove colour bar
    )

### Clouds
## Cloud isosurface
cloud_surf = go.Isosurface(
    z=Z.flatten(),
    x=X.flatten(),
    y=Y.flatten(),
    value=qc0_xyz.flatten(),
    isomin=.1,
    isomax=1,
    surface_count=4,
    opacity=0.6,
    showscale=False,
    )

### Theta
## theta on flat surface (model levels)
flat_theta_surf = go.Isosurface(
    z=Z[:,:,:topind_th].flatten(),
    x=X[:,:,:topind_th].flatten(),
    y=Y[:,:,:topind_th].flatten(),
    value=th0_xyz[:,:,:topind_th].flatten(),
    isomin=311,
    isomax=320,
    surface_count=4,
    opacity=0.7,
    colorscale='Hot',
    showscale=False
    )

## theta on topography (model levels + topog)
theta_surf = go.Isosurface(
    z=Zh[:,:,:topind_th].flatten(),
    x=X[:,:,:topind_th].flatten(),
    y=Y[:,:,:topind_th].flatten(),
    value=th0_xyz[:,:,:topind_th].flatten(),
    isomin=311,
    isomax=320,
    surface_count=4,
    opacity=0.7,
    colorscale='Hot',
    showscale=False
    )


### Wind
## wind quiver in 3d, tubes?
wind_tubes = go.Streamtube(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    u=u0_xyz.flatten(),
    v=v0_xyz.flatten(),
    w=w0_xyz.flatten(),
    colorscale='Blues',
    showscale=False,
    cmin=1,
    cmax=5, # for colorbar
    maxdisplayed=40, # how many tubes?
    sizeref=.5, #??
    )
## wind quiver almost using cones?
skipsize=20 # show 1 in skipsize cones along each dimension
skip=(slice(None,None,skipsize),slice(None,None,skipsize),slice(None,None,skipsize))
wind_cones = go.Cone(
    x=X[skip].flatten(),
    y=Y[skip].flatten(),
    z=Z[skip].flatten(),
    u=u0_xyz[skip].flatten(),
    v=v0_xyz[skip].flatten(),
    w=w0_xyz[skip].flatten(),
    colorscale='Blues',
    sizemode="absolute",
    sizeref=13,
    showscale=False,
    opacity=.7,
    #cmin=1,
    #cmax=5, # for colorbar
    #maxdisplayed=40, # how many tubes?
    )

## Surface level winds
sskip=(slice(None,None,skipsize),slice(None,None,skipsize),slice(None,None,None))
surf_wind_cones = go.Cone(
    x=X[sskip][:,:,0].flatten(),
    y=Y[sskip][:,:,0].flatten(),
    z=Z[sskip][:,:,15].flatten(),
    u=u0_xyz[sskip][:,:,0].flatten(),
    v=v0_xyz[sskip][:,:,0].flatten(),
    w=w0_xyz[sskip][:,:,0].flatten(),
    colorscale='Blues',
    sizemode="absolute",
    sizeref=2,
    showscale=False,
    opacity=.5,
    #cmin=1,
    #cmax=5, # for colorbar
    #maxdisplayed=40, # how many tubes?
    )

## Vertical winds contour at 500ish metres
updraft_height = 1700
upd_ind = np.sum(levh<updraft_height)
updrafts = go.Surface(
    z=Z[:,:,upd_ind],
    x=X[:,:,upd_ind],
    y=Y[:,:,upd_ind],
    colorscale='PiYG_r',
    surfacecolor=w0_xyz[:,:,upd_ind],
    opacity=.65,
    cmin=-2, 
    cmax=2,
    showscale=False, # remove colour bar,
    )

### SHOW FIGURES 
show_topog=False
show_cloud=False
show_theta=False
show_wind=False
show_combined=True

# fire surface
if show_topog:
    figf = go.Figure(data=topog_fire)
    figf.update_layout(northerly_view)
    #pio.show(figf)
    figf.write_image("test3d.png")
    
    
    
# cloud isosurface
if show_cloud:
    figqc = go.Figure(data=cloud_surf)
    figqc.update_layout(northerly_view)
    pio.show(figqc)

if show_theta:
    figth = go.Figure(data=theta_surf)
    figth.update_layout(northerly_view)
    pio.show(figth)

if show_wind:
    #figtube = go.Figure(data=wind_tubes)
    #pio.show(figtube)
    figcone = go.Figure(data=surf_wind_cones)
    pio.show(figcone)


if show_combined:
    # Combine theta and cloud isosurfaces, over topography surface
    figall = go.Figure(data=[flat_topog, updrafts,flat_theta_surf,cloud_surf])
    # Add trace of theta to longitude wall
    # figall.update_trace() #TODO
    figall.update_layout(standard_view)
    
    pio.show(figall)
    
    ## Try to save figure
    #figall.write_image('test3d.png')

# may need to restart orca...
pio.orca.shutdown_server()