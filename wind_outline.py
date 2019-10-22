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

def show_transects():
    """
    For sirivan and waroona show the transects on a contour map with final fire outline as well
    """
    
    # Show transects
    
    transects = plotting._transects_
    waroona_xticks = np.arange(115.8, 116.21, .1)
    waroona_yticks = np.arange(-33, -32.701, .05)
    sirivan_xticks = np.arange(149.2, 150.4, .2) #149.2, 150.4
    sirivan_yticks = np.arange(-32.4,-31.6,.1) #-32.4, -31.6
    for extentname, xticks, yticks in zip(['waroona','sirivan'],
                                          [waroona_xticks,sirivan_xticks],
                                          [waroona_yticks,sirivan_yticks]):
        extent = plotting._extents_[extentname]
        mr = '%s_run1'%extentname
        
        # read topog, fire
        topog = fio.read_topog(mr,extent=extent)
        #print("debug:",topog)
        lat,lon = topog.coord('latitude').points, topog.coord('longitude').points
        
        ff, = fio.read_fire(mr,extent=extent,
                            dtimes=[fio.model_outputs[mr]['filedates'][-1]])
        #print("DEBUG:",ff)
        plotting.map_topography(extent,topog.data,lat,lon)
        plt.title("Transects")
        plt.contour(lon,lat, np.transpose(ff[0].data),np.array([0]),colors='red')
        for transect in range(6):
            start,end = transects["%s%d"%(extentname,transect+1)]
            plt.plot([start[1],end[1]],[start[0],end[0], ], '--',
                     linewidth=2, label='X%d'%(transect+1),
                     marker='X', markersize=7,markerfacecolor='white')
        
        # add nearby towns
        plotting.map_add_locations_extent(extentname)
        
        plt.legend()
        plt.yticks(yticks)
        plt.xticks(xticks)
        pname="figures/transects_%s.png"%extentname
        fio.save_fig_to_path(pname,dpi=600,plt=plt)
    
    
def wind_outline(s,u,v,w,
                 qc, topog,
                 z,lat,lon,
                 transect,
                 ff=None,
                 extentname='waroona',
                 nquivers=18,
                 quiverscale=60,
                 ztop = 5000,
                 cloud_threshold=constants.cloud_threshold,
                ):
    '''
    311 Plot showing contour map and wind speed, along with near sites and transect
    312 plot showing vertical motion along transect
    313 plot showing wind speed along transect
    INPUTS:
        wind speed, x-wind, y-wind, vert-wind, cloud content, topography, model level altitude, lat,lon, datetime (for titles)
        ff: fire front array, negative where fire has burned, optional
        extentname = { 'waroona' | 'sirivan' } choice of two fire locations
        transect = [[lat0,lon0],[lat1,lon1]]
        nquivers: int, roughly how many quivers are desired
        quiverscale: int, changes how long the arrows are (also may need to fiddle)
        cloud_threshold: float, g/kg where cloud outline is drawn
    '''
    # set font sizes
    plotting.init_plots()
    
    # get plot extent, and transect
    extent = [lon[0],lon[-1],lat[0],lat[-1]]
    start,end = transect
    
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
             marker='>', markersize=7,markerfacecolor='white')
    
    # add nearby towns
    if extentname is not None:
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
                                                  ztop=ztop,
                                                  lines=None)
    plt.ylabel('height (m)')
    ## Add contour where clouds occur
    qcslice = utils.cross_section(qc,lat,lon,start,end, npoints=100)
    with warnings.catch_warnings():
        # ignore warning when there are no clouds:
        warnings.simplefilter('ignore')
        plt.contour(xslice,zslice,qcslice,np.array([cloud_threshold]),colors='teal')
    
    ax3 = plt.subplot(3,1,3)
    trs, trx, trz = plotting.transect_s(s,z,lat,lon,start,end,topog=topog, ztop=ztop)
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
    
    return plt

def outline_model_winds(model_run='sirivan_run1', hours=None, dpi=200):
    
    extentname=model_run.split('_')[0]
    extent=plotting._extents_[extentname]
    
    ## Read the cubes
    cubes = fio.read_model_run(model_run, fdtime=hours, extent=extent,
                               add_topog=True,
                               add_z=True, 
                               add_winds=True)
    
    zth, topog = cubes.extract(['z_th','surface_altitude'])
    u,v,s, qc, w = cubes.extract(['u','v','s','qc','upward_air_velocity'])
    cubetimes = utils.dates_from_iris(u)
    
    ## fire front
    ff, = fio.read_fire(model_run=model_run, dtimes=cubetimes, extent=extent, 
                        firefront=True)
    lat, lon = w.coord('latitude').points, w.coord('longitude').points
    
    ## loop over available timedate steps
    for i, dt in enumerate(cubetimes):
        stitle = dt.strftime("%Y %b %d %H:%M (UTC)")
        si, ui, vi, wi, qci, topogi, zthi, ffi = [s[i].data.data, u[i].data.data, 
                                                 v[i].data.data, w[i].data.data, 
                                                 qc[i].data.data, topog.data.data, 
                                                 zth.data.data, ff[i].data.data]

        ## loop over different transects
        for ii in range(6):
            transect = plotting._transects_["%s%d"%(extentname,ii+1)]
            plt = wind_outline(si, ui, vi, wi, qci, topogi, zthi, lat, lon,
                               transect=transect,
                               ff = ffi,
                               extentname=extentname)
            # Save figure into subfolder with transect identifier
            plt.suptitle(stitle)
            fio.save_fig(model_run, _sn_, cubetimes[i], plt,
                         subdir='transect_%d'%(ii+1),
                         dpi=dpi)

#########################################################################
#########################################################################
#########################################################################


if __name__ == '__main__':
    
    for mr in ['sirivan_run1']:#,'waroona_old','waroona_run1']:
        hours = fio.model_outputs[mr]['filedates']
        for hour in hours:
            print("info: wind_outline", mr, hour)
            outline_model_winds(mr, hours=[hour])
    