#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 15:14:27 2020
    Look at surfaces in 3d
        show surface level, and atmospheric stuff in a set of 3d axes
    may require extra modules to be loaded...
    module load gtk
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
# Make browser the default plotly renderer
pio.renderers.default = "browser"

## Run these if running on local laptop:
import sys
if "g/data" not in sys.prefix:
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
_verbose_ = True

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

###
## METHODS
###

def verbose(*args):
    if _verbose_:
        print("INFO(verbose):",*args)

def cube_to_xyz(cube,
                ztopind=-1):
    """
    take iris cube [lev, lat, lon]
    pull out the data and reshape it to [lon,lat,lev]
    """
    assert len(cube.shape)==3, "cube is not 3-D"
    data = cube[:ztopind,:,:].data.data
    # data is now a level, lat, lon array... change to lon, lat, lev
    xyz = np.moveaxis(np.moveaxis(data,0,2),0,1)
    return xyz
    

def create_figure(gofigures, 
                  camera_eye=[1.4,-1.6,.35], 
                  aspectratio=[1,1,0.8], 
                  filename=None,
                  **layoutargs):
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
                                 camera_eye=dict(x=cx,y=cy,z=cz),
                                 xaxis = dict(title='lon'),
                                 yaxis = dict(title='lat'),
                                 zaxis = dict(title='alt'),
                                 ),
                      #margin=dict(#l=0,
                      #            #r=0,
                      #            t=0.5,
                      #            b=0.2,
                      #            ),
                      )
    
    if len(layoutargs) > 0:
        fig.update_layout(**layoutargs)
    
    if filename is None:
        pio.show(fig)
    else:
        ## Try to save figure
        fio.make_folder(filename)
        print("INFO: Saving Image: ",filename)
        fig.write_image(filename)

def cloud_system(model_run='waroona_run2', hour=20, 
                theta_height=1000, theta_min=311, theta_max=320,
                vert_motion_height = 1700,
                top_height=8000, send_to_browser=False,
                extent=None,
                HSkip=5):
    """
    Read an hour of model output, plot it in 3d using plotly
    saves output as .png
    ARGUMENTS:
        hour: which output hour (first is 0, last is 23 or -1)
        theta_height: how high are we looking regarding potential temp?
        theta_min, theta_max: min and max potential temperature to draw isosurface
        vert_motion_height: altitude of vertical motion surface,
        height_top: how high (m) included in plot
        send_to_browser: instead of trying to save figures, send one to the browser (interactive)
    """
    
    hours=fio.model_outputs[model_run]['filedates']
    dtime=hours[hour]
    if extent is None:
        extentname = model_run.split('_')[0]
        extent = plotting._extents_[extentname]
    
    cubes = fio.read_model_run(model_run, 
                               fdtime=dtime, 
                               extent=extent, 
                               add_theta=True, 
                               add_topog=True, 
                               add_winds=True,
                               HSkip=HSkip)
    
    th, qc = cubes.extract(['potential_temperature','qc'])
    # datetimes in hour output
    cubetimes = utils.dates_from_iris(th)
    
    ff, = fio.read_fire(model_run,
                        dtimes=cubetimes,
                        extent=extent,
                        HSkip=HSkip)
    
    # Get the rest of the desired data
    topog, = cubes.extract(['surface_altitude'])
    d_topog = topog.data.data.T # convert to lon,lat
    u,v,w = cubes.extract(['u','v','upward_air_velocity'])
    levh  = qc.coord('level_height').points
    topind = np.sum(levh<top_height)
    topind_th = np.sum(levh<theta_height)
    # these are level, lat, lon cubes
    lat,lon = qc.coord('latitude').points, qc.coord('longitude').points
    
    # dimensional mesh
    X,Y,Z = np.meshgrid(lon,lat,levh) 
    ## X Y Z are now [lat, lon, lev] for some reason
    [X,Y,Z] = [ np.moveaxis(arr,0,1) for arr in [X,Y,Z]]
    ## Now they are lon, lat, lev
    ## Cut down to desired level
    [X, Y, Z] = [ arr[:,:,:topind] for arr in [X,Y,Z]]

    # topography surface
    topog_layer = go.Surface(
        z=Z[:,:,0],
        x=X[:,:,0],
        y=Y[:,:,0],
        colorscale='earth', # was not reversed on local laptop
        reversescale=True,
        surfacecolor=d_topog,
        showscale=False, # remove colour bar,
    )
    
    namedlocs=[]
    namedlocs_lats = []
    namedlocs_lons = []
    for (namedloc, (loclat, loclon)) in plotting._latlons_.items():
        #print(namedloc, loclat, loclon)
        if loclon < extent[1] and loclon > extent[0] and loclat < extent[3] and loclat > extent[2]:
            if 'fire' not in namedloc and 'pyrocb' not in namedloc:
                namedlocs.append(namedloc)
                namedlocs_lats.append(loclat)
                namedlocs_lons.append(loclon)
    
    for hi, cubetime in enumerate(cubetimes):
        verbose("Creating surfaces")
        
        # get cloud, theta, vert motion, firefront in terms of lon,lat,lev
        d_qc = cube_to_xyz(qc[hi],ztopind=topind)
        d_th = cube_to_xyz(th[hi],ztopind=topind_th)
        d_w = cube_to_xyz(w[hi],ztopind=topind)
        
        #d_ff = ff[hi].data.data # firefront already in lon,lat shape
        
        # surfaces to be plotted in 3d
        surface_list = [topog_layer]
        
        ## Points for waroona, yarloop
        locations_scatter = go.Scatter3d(
            x=namedlocs_lons,
            y=namedlocs_lats,
            z=[0]*len(namedlocs),
            mode='markers',
            marker=dict(
                size=12,
                color='black',           # array/list of desired values
                #colorscale='Viridis',   # choose a colorscale
                opacity=0.8
                ),
            )
        surface_list.append(locations_scatter)
        
        ## atmospheric heat (theta)
        # TODO: smaller markers, labelled by location name
        if np.sum(d_th > theta_min) > 0:
            theta_surf = go.Isosurface(
                z=Z[:,:,:topind_th].flatten(),
                x=X[:,:,:topind_th].flatten(),
                y=Y[:,:,:topind_th].flatten(),
                value=d_th.flatten(),
                isomin=theta_min,
                isomax=theta_max,
                surface_count=4,
                opacity=0.7,
                colorscale='Hot',
                showscale=False
                )
            surface_list.append(theta_surf)
            verbose("adding heat surface")
        
        ## Vertical winds contour at some altitude 
        vm_ind = np.sum(levh<vert_motion_height)
        vert_motion_layer = go.Surface(
            z=Z[:,:,vm_ind],
            x=X[:,:,vm_ind],
            y=Y[:,:,vm_ind],
            colorscale='PiYG', # This was PiYG_r on local laptop
            reversescale=True,  # should be equivalent to _r
            surfacecolor=d_w[:,:,vm_ind],
            opacity=.65,
            cmin=-2, 
            cmax=2,
            showscale=False, # remove colour bar,
            )
        surface_list.append(vert_motion_layer)
        
        ## Cloud isosurface
        cloud_surf = go.Isosurface(
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
        surface_list.append(cloud_surf)
        
        #
        ### 3d figure:
        #
        
        ## title, lables
        layoutargs = dict(
            title=cubetime.strftime('Clouds %Y%m%d%H%M(UTC)'),
            font=dict(
                    family="Courier New, monospace",
                    size=18,
                    color="#222222"
                ),
            )

        
        figname = None
        if not send_to_browser:
            #figname = cubetime.strftime('figures/threedee/test_%Y%m%d%H%M.png')
            figname = fio.standard_fig_name(model_run,_sn_,cubetime)
        create_figure(surface_list, filename=figname, **layoutargs)

def downslope_system(model_run = 'waroona_run3', 
                     hour=20,
                     extent=[115.88, 116.0, -32.86,-32.83],
                     HSkip=None,
                     send_to_browser=True):
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

if __name__=='__main__':
    wider_waroona = plotting._extents_['waroona']
    wider_waroona[0] -= .4 # move west edge west
    wider_waroona[1] += .2 # move east edge east
    wider_waroona[2] -= .1 # move south edge south
    
    # check downslope stuff!
    if True:
        overfire=[116.09, 116.19, -32.91,-32.84]
        downslope_system(model_run='waroona_run1',hour=18,
                         extent=overfire,
                         HSkip=None)
    
    
    if False:
        
        # Save a bunch of images
        for hour in [17]:#range(15,24):
            
            #theta_min=311
            #theta_max=320
            #vert_motion_height = 1700
            cloud_system(model_run='waroona_run2',
                        hour = hour,
                        extent=wider_waroona,
                        HSkip=8,
                        top_height=13500,
                        theta_height=2000,
                        send_to_browser=True,)
