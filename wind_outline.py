#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 20:45:05 2019
    Winds outline script
@author: jesse
"""

import matplotlib
# don't plot on screen, send straight to file
# this is for lack of NCI display
matplotlib.use('Agg',warn=False)

# plotting stuff
import matplotlib.colors as col
import matplotlib.pyplot as plt
#import matplotlib.ticker as tick
#import matplotlib.patches as mpatches
import numpy as np
from datetime import datetime,timedelta
import warnings

# local modules
from utilities import plotting, utils, fio, constants

###
## GLOBALS
###
_sn_ = 'wind_outline'


def wind_outline(s,u,v,w,
                 qc, topog,
                 z,lat,lon,
                 dtime,
                 ff=None,
                 extentname='waroona',
                 transect=1, 
                 nquivers=18,
                 quiverscale=60,
                 cloud_threshold=constants.cloud_threshold,
                ):
    '''
    311 Plot showing contour map and wind speed, along with near sites and transect
    312 plot showing vertical motion along transect
    313 plot showing wind speed along transect
    INPUTS:
        dictionary with topog, vert motion, wind speed, wind direction, zth,lat,lon, time
        timestep: data will have a time dimension
        extentname = { 'waroona' | 'sirivan' } choice of two fire locations
        transect = int from 1 to 6 for transect choice
        vectorskip reduces quiver density (may need to fiddle)
        quiverscale changes how long the arrows are (also may need to fiddle)
        ext is the plot extension { '.png' | '.eps' }
    '''
    stitle = dtime.strftime("%Y %b %d %H:%M (UTC)")
    
    # set font sizes
    plotting.init_plots()
    
    # get plot extent, and transect
    extent = plotting._extents_[extentname]
    start,end = plotting._transects_["%s%d"%(extentname,transect)]
    
    plt.figure(figsize=[7,10])
    plt.subplot(3,1,1)
    
    # top panel is topography
    plotting.map_topography(extent,topog,lat,lon)
    plt.title('Topography, winds')
    
    # Add fire front contour
    if ff is not None:
        with warnings.catch_warnings():
            # ignore warning when there are no fires:
            warnings.simplefilter('ignore')
            plt.contour(lon,lat,np.transpose(ff),np.array([0]), 
                        colors='red')
    
    # start to end x=[lon0,lon1], y=[lat0, lat1]
    plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
             linewidth=2, 
             marker='X', markersize=7,markerfacecolor='white')
    
    # add nearby towns
    plotting.map_add_locations_extent(extentname)
    
    # Quiver, reduce arrow density
    vsu = len(lon)//nquivers
    vsv = len(lat)//nquivers
    skip = (slice(None,None,vsv),slice(None,None,vsu))
    
    ## colour the arrows
    # map wind speed to colour map domain [0, 1]
    norm = col.Normalize()
    norm.autoscale(s[skip])
    plt.get_cmap(plotting._cmaps_['windspeed'])
    plt.quiver(lon[skip[1]],lat[skip[0]],u[0][skip],v[0][skip], scale=quiverscale)
               #color=cmap(norm(s[skip])), 
    
    ## Second row is transect plots
    plt.subplot(3,1,2)
    wslice, xslice, zslice  = plotting.transect_w(w,z, lat, lon,start,end,
                                                  npoints=100,topog=topog,
                                                  lines=None)
    plt.ylabel('height (m)')
    ## Add contour where clouds occur
    qcslice = utils.cross_section(qc,lat,lon,start,end, npoints=100)
    with warnings.catch_warnings():
        # ignore warning when there are no clouds:
        warnings.simplefilter('ignore')
        plt.contour(xslice,zslice,qcslice,np.array([cloud_threshold]),colors='teal')
    
    ax3 = plt.subplot(3,1,3)
    trs, trx, trz = plotting.transect_s(s,z,lat,lon,start,end,topog=topog)
    xticks,xlabels = plotting.transect_ticks_labels(start,end)
    plt.xticks(xticks[0::2],xlabels[0::2]) # show transect start and end
    
    # Annotate max wind speed
    # only care about lower troposphere
    # 60 levels is about 2800m, 80 levels is about 4700m, 70 levels : 3700m
    upto = 70 
    mloc = np.unravel_index(np.argmax(trs[:upto,:],axis=None),trs[:upto,:].shape)
    note="max = %5.1f m/s\n  (at %4.0f m) "%(trs[:upto,:][mloc], trz[:upto,:][mloc])
    trans = ax3.get_xaxis_transform() # x in data untis, y in axes fraction
    ax3.annotate(note, fontsize=15,
                 xy=(0.33, -0.2 ), xycoords=trans)
    
    # Save figure into animation folder with numeric identifier
    plt.suptitle(stitle)
    return plt


def winds_outline_hour(dtime, model_version='waroona_run1', dpi=300):
    '''
    Loop over transects over waroona and make the wind outline figures
    '''
    extentname=model_version.split('_')[0]
    extent=plotting._extents_[extentname]
    
    # Read the cubes
    cubes = fio.read_model_run(model_version, fdtime=dtime, extent=extent,
                               add_topog=True,
                               add_z=True, 
                               add_winds=True)
    
    zth, topog = cubes.extract(['z_th','surface_altitude'])
    u,v,s, qc, w = cubes.extract(['u','v','s','qc','upward_air_velocity'])

    cubetimes = utils.dates_from_iris(u)
    
    ## fire front
    ff1, = fio.read_fire(dtimes=cubetimes, extent=extent, firefront=True)
    lat,lon = w.coord('latitude').points, w.coord('longitude').points
    
    # also loop over different transects
    for i in range(6):
        for tstep in range(len(cubetimes)):
            ff = None
            if model_version == 'waroona_run1':
                ff = ff1[tstep].data.data
            plt = wind_outline(s[tstep].data.data, u[tstep].data.data, v[tstep].data.data, w[tstep].data.data,
                               qc[tstep].data.data, topog.data.data,
                               zth.data.data, lat,lon,
                               dtime=cubetimes[tstep],
                               ff = ff,
                               extentname=extentname,
                               transect=i+1)
            fio.save_fig(model_version,_sn_,cubetimes[tstep],plt, 
                         subdir='transect_%d'%i,
                         dpi=dpi)

#########################################################################
#########################################################################
#########################################################################


if __name__ == '__main__':
    
    for dtime in [ datetime(2016,1,6,7) + timedelta(hours=x) for x in range(2) ]:
        for mv in ['waroona_old','waroona_run1']:
            winds_outline_hour(dtime, mv)
    