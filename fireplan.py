#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 16:28:18 2019
    Show fire spread and intensity over time
@author: jesse
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
from datetime import datetime, timedelta


from utilities import plotting, utils, fio 

###
## GLOBALS
###
_sn_ = 'fireplan'
PFT = {'waroona_run1':{'data':np.array([34.6]), # manually calculated PFT
                       'units':'Gigawatts',
                       'time':np.array([datetime(2016,1,6,6)]), # calculated at these times in UTC
                       'latlon':[-32.89 -0.004, 116.17+0.009], # ~ 1km from fire
                       'latlon_stamp':'fire_waroona_upwind',
                       },
       'waroona_old':{'data':np.array([138.7]), # manually calculated PFT
                      'units':'Gigawatts',
                      'time':np.array([datetime(2016,1,6,6)]), # calculated at these times in UTC
                      'latlon':[-32.89 -0.004, 116.17+0.009], # ~ 1km from fire
                      'latlon_stamp':'fire_waroona_upwind',
                       },
      }


def fireplan(ff, extentname='waroona', fire_contour_map = 'autumn',
             show_cbar=True, cbar_XYWH= [0.2,0.6,.3,.02],
             fig=None,subplotxyn=None,draw_gridlines=False,gridlines=None):
    '''
    show google map of extent, with fire front over time overplotted
    
    ARGUMENTS:
        ff: iris.cube.Cube with time, longitude, latitude dimensions
        extentname: 'waroona' or 'sirivan'
        fire_contour_map: how will firefront contours be coloured
        show_cbar: bool
            draw a little colour bar showing date range of contours?
        cbar_XYHW: 4-length list where to put cbar
        fig,... arguments to plotting.map_google()
    '''
    lon,lat = ff.coord('longitude').points, ff.coord('latitude').points
    nt,nx,ny = ff.shape
    extent = plotting._extents_[extentname] 
    
    # Get datetimes from firefront cube
    fire_front_hours = utils.dates_from_iris(ff) 
    # How many fire contours do we have?
    minff = np.min(ff.data,axis=(1,2))
    subset_with_fire = minff < 0
    ff_f = ff[subset_with_fire]
    fire_front_hours = fire_front_hours[subset_with_fire]
    nt = np.sum(subset_with_fire)
    
    ## fire contour colour map
    cmap = matplotlib.cm.get_cmap(fire_contour_map)
    # we have 14 hours with fire
    rgba = cmap(np.linspace(0,1,nt))
    
    ## PLOTTING
    
    # First plot google map 
    plotting.map_google(extent, zoom=11, fig=fig, subplotxyn=subplotxyn, 
                        draw_gridlines=draw_gridlines, gridlines=gridlines)
    
    #plotting.map_add_locations(['waroona'],text=['Waroona'],dx=.03,dy=-.006)
    #plotting.map_add_locations(['fire_waroona'],text=['F160'],dx=-.001,dy=-.006)
    
    # plot contours and label time stamp
    utcstamp=[]
    for ii,dt in enumerate(fire_front_hours):
        
        utcstamp.append(dt.strftime('%b %d, %HH(UTC)'))
        
        ffdata = ff_f[ii].data.data
        assert np.min(ffdata) < 0, 'no fire yet'    
        cs = plt.contour(lon,lat,ffdata.T, np.array([0]), 
                         colors=[rgba[ii]],linewidths=1)
        
        ## Add dstamp as inline label...
        #if dt.hour%6 == 0:
        #    plt.clabel(cs, inline=True, fmt=dstamp, fontsize=9, fontcolor='k')
        
    ## Add tiny colour bar showing overall time of fire
    if show_cbar:
        cax = fig.add_axes(cbar_XYWH)
        #cax.axes.get_xaxis().set_visible(False)
        cax.axes.get_yaxis().set_visible(False)
        cbar = matplotlib.colorbar.ColorbarBase(cax,cmap=plt.get_cmap(fire_contour_map), orientation='horizontal')
        cbar.set_ticks([0,1])
        cbar.set_ticklabels([utcstamp[0],utcstamp[-1]])
        
    
def fireplan_waroona():
    '''
    Show fire outline over time
    
    ARGUMENTS:
        extentname: which extent are we looking at
        front_start: from what hour do we have a fire front?
    '''
    # Read fire output
    extentname = 'waroona'
    front_start=datetime(2016,1,6,1)
    extent = plotting._extents_[extentname+'z'] # Super zoomed
    
    fire_front_hours = [front_start + timedelta(hours=x) for x in range(15)]
    FFront, SHeat, FSpeed = fio.read_fire(dtimes=fire_front_hours, extent=extent,
                                          firefront=True, sensibleheat=True, firespeed=True)
    lon,lat = FFront.coord('longitude').points, FFront.coord('latitude').points
    
    ## fire contour colour map
    fire_contour_map = 'autumn'
    cmap = matplotlib.cm.get_cmap(fire_contour_map)
    # we have 14 hours with fire
    rgba = cmap(np.linspace(0,1,len(fire_front_hours)))
    
    ## PLOTTING
    fig = plt.figure(figsize=[8,9])
    ax1 = plt.subplot(2,1,1)
    
    # First plot topography
    topog = fio.read_topog('waroona_run1',extent=extent)
    plt.contourf(lon,lat, topog.data, 100, cmap='terrain', vmin=-100)
    plotting.map_add_locations(['waroona'],text=['Waroona'],dx=.03,dy=-.006)
    plotting.map_add_locations(['fire_waroona_upwind'],text=['F160'],dx=-.001,dy=-.006)
    
    # plot contours and label time stamp
    utcstamp=[]
    for ii,dt in enumerate(fire_front_hours):
        
        utcstamp.append(dt.strftime('%b %d, %HH(UTC)'))
        #dstamp = (dt+LT_offset).strftime("%H(LT)")
        FF = FFront[ii].data.data
        assert np.min(FF) < 0, 'no fire yet'    
        cs = plt.contour(lon,lat,FFront[ii].data.T, np.array([0]), 
                         colors=[rgba[ii]],linewidths=1)
        
        ## Add dstamp as inline label...
        #if dt.hour%6 == 0:
        #    plt.clabel(cs, inline=True, fmt=dstamp, fontsize=9, fontcolor='k')
        
    ## Add tiny colour bar showing overall time of fire   
    cax = fig.add_axes([0.2,0.6,.3,.02])
    #cax.axes.get_xaxis().set_visible(False)
    cax.axes.get_yaxis().set_visible(False)
    cbar = matplotlib.colorbar.ColorbarBase(cax,cmap=plt.get_cmap(fire_contour_map), orientation='horizontal')
    cbar.set_ticks([0,1])
    cbar.set_ticklabels([utcstamp[0],utcstamp[-1]])
    
    ## subplot 2
    ax2 = plt.subplot(4,1,3)
    maxflux=np.max(SHeat.data,axis=0)
    cs=plt.contourf(lon, lat, maxflux.T, 30, 
                    cmap='inferno', 
                    locator=ticker.LogLocator())
    plt.colorbar(cs, orientation='horizontal', pad=0)
    plt.title('sensible heat flux (?)',y=.74)
    plt.xticks([],[])
    plt.yticks([],[])
    
    ax3 = plt.subplot(4,1,4)
    maxspeed=np.max(FSpeed.data,axis=0)
    cs=plt.contourf(lon, lat, maxspeed.T, 30, 
                    cmap=plotting._cmaps_['windspeed'])
    plt.colorbar(cs, orientation='horizontal', pad=0)
    plt.title('Max Firespeed (m/s?)', y=.74)
    plt.xticks([],[])
    plt.yticks([],[])
    
    fio.save_fig('waroona_run1',_sn_, 'fire_spread', plt)