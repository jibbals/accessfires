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

## Run these if running on local laptop:
import sys
if "g/data" not in sys.prefix:
    # Make browser the default plotly renderer
    pio.renderers.default = "browser"
    # turn on orca (server that takes plotly interactive images and saves them to static png)
    pio.orca.ensure_server()
    # check orca status
    pio.orca.status

# local modules
from utilities import plotting, utils, fio

###
## GLOBALS
###
_sn_ = 'threedee'


###
## METHODS
###

def cube_to_xyz(cube,
                ztop=-1):
    """
    take iris cube [lev, lat, lon]
    pull out the data and reshape it to [lon,lat,lev]
    """
    assert len(cube.shape)==3, "cube is not 3-D"
    data = cube[:ztop,:,:].data.data
    # data is now a level, lat, lon array... change to lon, lat, lev
    xyz = np.moveaxis(np.moveaxis(data,0,2),0,1)
    return xyz
    

def create_figure(gofigures, 
                  camera_eye=[1.5,-1.5,.8], 
                  aspectratio=[1,1,1], 
                  filename=None):
    """
    show or save a figure
    gofigures: list of go objects (e.g., [go.Surface(...),go.Isosurface(...)])
    view: where the camera viewpoint is located
    aspectratio: the aspect ratio (todo: figure out)
    filename: if not None, try to save figure 
    """
    # make figure
    fig = go.Figure(data=gofigures)
    # place camera
    cx,cy,cz=camera_eye
    ax,ay,az=aspectratio
    fig.update_layout(scene=dict(aspectratio=dict(x=ax,y=ay,z=az),
                                 camera_eye=dict(x=cx,y=cy,z=cz)))
    if filename is None:
        pio.show(fig)
    else:
        ## Try to save figure
        fig.write_image(filename)

def save_system(model_run='waroona_run2', hour=20, 
                theta_height=1000, theta_min=311, theta_max=320,
                top_height=8000, send_to_browser=False):
    """
    Read an hour of model output, plot it in 3d using plotly
    saves output as .png
    ARGUMENTS:
        hour: which output hour (first is 0, last is 23 or -1)
        theta_height: how high are we looking regarding potential temp?
        theta_min, theta_max: min and max potential temperature to draw isosurface
        height_top: how high (m) included in plot
        send_to_browser: instead of trying to save figures, send one to the browser (interactive)
    """
    
    hours=fio.model_outputs[model_run]['filedates']
    dtime=hours[hour]
    extentname = model_run.split('_')[0]
    extent = plotting._extents_[extentname]
    
    cubes = fio.read_model_run(model_run, fdtime=dtime, extent=extent, 
                               add_theta=True, add_topog=True, add_winds=True)
    ff, = fio.read_fire(model_run,dtimes=[dtime],extent=extent)
    d_ff = ff[0].data.data # firefront already in lon,lat shape
    
    th, qc = cubes.extract(['potential_temperature','qc'])
    topog, = cubes.extract(['surface_altitude'])
    d_topog = topog.data.data.T # convert to lon,lat
    u,v,w = cubes.extract(['u','v','upward_air_velocity'])
    levh  = qc.coord('level_height').points
    topind = np.sum(levh<top_height)
    topind_th = np.sum(levh<theta_height)
    # these are level, lat, lon cubes
    lat,lon = qc.coord('latitude').points, qc.coord('longitude').points
    cubetimes = utils.dates_from_iris(th)
    
    # dimensional mesh
    X,Y,Z = np.meshgrid(lon,lat,levh) 
    ## X Y Z are now [lat, lon, lev] for some reason
    [X,Y,Z] = [ np.moveaxis(arr,0,1) for arr in [X,Y,Z]]
    ## Now they are lon, lat, lev
    
    # topography surface
    surface_topog = go.Surface(
        z=Z[:,:,0],
        x=X[:,:,0],
        y=Y[:,:,0],
        colorscale='earth',
        surfacecolor=d_topog,
        showscale=False, # remove colour bar,
    )
    
    # get cloud, theta, vert motion, firefront in terms of lon,lat,lev
    d_qc = cube_to_xyz(qc[0],ztop=topind)
    d_th = cube_to_xyz(th[0],ztop=topind_th)
    d_w = cube_to_xyz(w[0],ztop=topind)
    d_levh = levh[:topind]
    
    ## Cloud isosurface
    isosurface_cloud = go.Isosurface(
        z=Z.flatten(),
        x=X.flatten(),
        y=Y.flatten(),
        value=d_qc.flatten(),
        isomin=.1,
        isomax=1,
        surface_count=3,
        opacity=0.6,
        showscale=False,
    )
    
    #
    ### 3d figure:
    #
    figname = None
    if not send_to_browser:
        figname = time.strftime('figures/threedee/test_%Y%m%d%H%M.png')
    create_figure([surface_topog,isosurface_cloud], filename=figname)

if __name__=='__main__':
    print("TESTING THREEDEE")
    save_system(send_to_browser=True)
