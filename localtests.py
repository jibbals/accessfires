#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""


# pretend we are in threedee module (this code will go there when ready)
from threedee import *
_sn_ = 'localtests'

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

overfire=[116.09, 116.19, -32.91,-32.84]
model_run='waroona_run1'
hour=18
extent=overfire
extent=[115.88, 116.0, -32.86,-32.83]
HSkip=2
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
u,v,w = cubes.extract(['u','v','upward_air_velocity'])
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

Z_3d = np.zeros(np.shape(Z))
for i in range(Z.shape[2]): # lon,lat,lev - loop over levels
    Z_3d[:,:,i] = Z[:,:,i] + d_topog # add topography to each level

for ti, cubetime in enumerate(cubetimes):
    surface_list=[]
    
    # surface sensible heat flux [t,lon,lat]
    d_sh = sh[ti].data.data
    
    d_zth = z_th[:topind,:,:].data.data # [lev, lat, lon]
    d_th = th[ti,:topind,:,:].data.data # [lev, lat, lon]
    #d_th = threedee.cube_to_xyz(th[ti],ztopind=topind) # [lon,lat,lev]
    
    d_w = cube_to_xyz(w[ti],ztopind=topind) # lon, lat, lev
    
    # topography surface
    topog_layer = go.Surface(
        z=Z[:,:,0],
        x=X[:,:,0],
        y=Y[:,:,0],
        colorscale='earth', # surface heat colormap 
        #cmin=shcmin,
        #cmax=shcmax,
        reversescale=True,
        surfacecolor=d_topog, # colour by sensible heat
        #opacityscale=[[0.0, 0], [100.0, .8], [shcmax, 1]], 
        #hidesurface=True,
        showscale=False, # remove colour bar,
        #coloraxis='coloraxis', # first colour bar
    )
    surface_list.append(topog_layer)
    
    ## volume plot showing vertical motion 
    wmax = 5
    wfade= 1
    w_volume=go.Volume(
        x=X.flatten(), 
        y=Y.flatten(), 
        z=Z.flatten(),
        value=d_w.flatten(),
        isomin=-1*wmax,
        isomax=wmax,
        opacity=0.1, # max opacity
        opacityscale=[[-1*wmax, 1], [-1*wfade, 0], [wfade, 0], [wmax, 1]], # only show |w|>wfade
        surface_count=20,
    )
    surface_list.append(w_volume)
    ### Paint figure
    figname = None
    if not send_to_browser:
        #figname = cubetime.strftime('figures/threedee/test_%Y%m%d%H%M.png')
        figname = fio.standard_fig_name(model_run,_sn_,cubetime,subdir='downslope')
    
    layoutargs = dict(title=cubetime.strftime('%Y%m%d%H%M(UTC)'),
                      font=dict(size=18, color="#111111"),)
                                
    create_figure(surface_list, filename=figname, **layoutargs)


