#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 16:28:18 2019
    Show fire spread and intensity over time
@author: jesse
"""

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import ticker, colors, patches
import numpy as np
from datetime import datetime, timedelta
import iris

from cartopy import crs as ccrs

from utilities import plotting, utils, fio 

###
## GLOBALS
###
_sn_ = 'fireplan'
PFT = {'waroona_run1':{'data':np.array([42.1,34.6]), # manually calculated PFT
                       'units':'Gigawatts',
                       'time':np.array([datetime(2016,1,6,5,10),datetime(2016,1,6,6)]), # calculated at these times in UTC
                       'latlon':[-32.89 -0.004, 116.17+0.009], # ~ 1km from fire
                       'latlon_stamp':'fire_waroona_upwind',
                       },
       'waroona_old':{'data':np.array([78.4,138.7]), # manually calculated PFT
                      'units':'Gigawatts',
                      'time':np.array([datetime(2016,1,6,5),datetime(2016,1,6,6)]), # calculated at these times in UTC
                      'latlon':[-32.89 -0.004, 116.17+0.009], # ~ 1km from fire
                      'latlon_stamp':'fire_waroona_upwind',
                       },
      }


def fireplan(ff, extentname='waroonaz', fire_contour_map = 'autumn',
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
    crs_data = ccrs.PlateCarree()
    
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
    rgba = cmap(np.linspace(0,1,nt))
    
    ## PLOTTING
    
    # First plot google map 
    gfig, gax, gproj = plotting.map_google(extent, zoom=12, 
                                           fig=fig, subplotxyn=subplotxyn, 
                                           draw_gridlines=draw_gridlines, 
                                           gridlines=gridlines)
    
    #plotting.map_add_locations(['waroona'],text=['Waroona'],dx=.03,dy=-.006)
    #plotting.map_add_locations(['fire_waroona'],text=['F160'],dx=-.001,dy=-.006)
    
    # plot contours and label time stamp
    utcstamp=[]
    for ii,dt in enumerate(fire_front_hours):
        
        utcstamp.append(dt.strftime('%b %d, %HH(UTC)'))
        
        ffdata = ff_f[ii].data.data
        assert np.min(ffdata) < 0, 'no fire yet'    
        plt.contour(lon,lat,ffdata.T, np.array([0]), 
                    colors=[rgba[ii]],linewidths=1, 
                    transform=crs_data)
        
    ## Add tiny colour bar showing overall time of fire
    if show_cbar:
        cax = fig.add_axes(cbar_XYWH)
        #cax.axes.get_xaxis().set_visible(False)
        cax.axes.get_yaxis().set_visible(False)
        cbar = matplotlib.colorbar.ColorbarBase(cax,cmap=plt.get_cmap(fire_contour_map), orientation='horizontal')
        cbar.set_ticks([0,1])
        cbar.set_ticklabels([utcstamp[0],utcstamp[-1]])
        plt.sca(gax) # return focus to newly created plot
    
    return gfig, gax, gproj

def fireplan_waroona():
    '''
    Show fire outline over time at waroona
    '''
    # Read fire output
    extentname = 'waroonaz' # fire zoomed
    #front_start=datetime(2016,1,6,1)
    extent = plotting._extents_[extentname] 
    
    # just read hourly
    first_day = [datetime(2016,1,5,15) + timedelta(hours=x) for x in range(24)]
    FFront, SHeat, FSpeed = fio.read_fire(dtimes=first_day, extent=extent,
                                          firefront=True, sensibleheat=True, firespeed=True)
    lon,lat = FFront.coord('longitude').points, FFront.coord('latitude').points
        
    ## PLOTTING
    
    ## First plot fire front contours over google map in a 211 subplot
    fig = plt.figure(figsize=[8,9])
    
    fireplan(FFront, extentname=extentname,
             fig=fig, subplotxyn=[2,1,1],
             show_cbar=True, cbar_XYWH=[.2,.6,.3,.02])
    
    plotting.map_add_locations(['fire_waroona_upwind'],text=['F160'],dx=-.001,dy=-.006, proj=ccrs.PlateCarree())
    
    ## subplot 2
    ax2 = plt.subplot(4,1,3)
    maxflux=np.sum(SHeat.data.data,axis=0) + 0.01 # get rid of zeros
    levels = np.sort(np.union1d(np.power(10,np.arange(2,6)),5*np.power(10,np.arange(2,6))))
    cs=plt.contourf(lon, lat, maxflux.T,
                    levels, # color levels I think...
                    norm=colors.LogNorm(),
                    vmin=100,
                    cmap='gnuplot2_r',
                    #locator=ticker.LogLocator(),
                    )
    plt.colorbar(cs, orientation='horizontal', pad=0, extend='max')
    plt.title('Total sensible heat flux (W/m2 ?)',y=.74)
    plt.xticks([],[])
    plt.yticks([],[])
    
    ax3 = plt.subplot(4,1,4)
    maxspeed=np.max(FSpeed.data,axis=0)
    cs=plt.contourf(lon, lat, maxspeed.T, 30, 
                    cmap=plotting._cmaps_['windspeed'])
    plt.colorbar(cs, orientation='horizontal', pad=0)
    plt.title('Max Firespeed (m/s ?)', y=.74)
    plt.xticks([],[])
    plt.yticks([],[])
    
    fio.save_fig('waroona_run1',_sn_, 'fire_spread', plt)
    
def fire_power_waroona():
    
    # Read fire output
    extentname = 'waroonaz' # fire zoomed
    #front_start=datetime(2016,1,6,1)
    extent = plotting._extents_[extentname] 
    
    # read all the fire data
    ff, sh = fio.read_fire(dtimes=None, extent=extent, firefront=True, sensibleheat=True)
    lon,lat = ff.coord('longitude'), ff.coord('latitude')
    # just read hourly for fireplan
    ftimes = utils.dates_from_iris(ff)
    hourinds = [ft.minute==0 for ft in ftimes]
    hours = ftimes[hourinds]
    
    ## First plot fire front contours over google map in a 211 subplot
    fig = plt.figure(figsize=[8,8])
    
    _, ax, gproj = fireplan(ff[hourinds], extentname=extentname,
                            fig=fig, subplotxyn=[2,1,1],
                            show_cbar=True, cbar_XYWH=[.2,.6,.3,.02])
    
    plotting.map_add_locations(['fire_waroona_upwind'],text=['F160'],dx=-.001,dy=-.006, proj=ccrs.PlateCarree())
    plt.title('Hourly fire front contour')
    ## ADD SUBSET RECTANGLE
    #    ax.add_patch(patches.Rectangle(xy=botleft,
    #                                   width=width, 
    #                                   height=width,
    #                                   #facecolor=None,
    #                                   fill=False,
    #                                   edgecolor='red',
    #                                   linewidth=2,
    #                                   #linestyle='-',
    #                                   alpha=0.9,
    #                                   transform=cartopy.crs.PlateCarree()
    #                                   ))
    
    ### get areas in m2
    # Add boundaries to grid
    lat.guess_bounds()
    lon.guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(ff)
    
    firepower = sh.data.data * grid_areas # W/m2 * m2
    firepower = firepower/1e9 # Watts to Gigawatts
    # also PFT values:
    run1_xy = PFT['waroona_run1']['time'],PFT['waroona_run1']['data']
    old_xy = PFT['waroona_old']['time'],PFT['waroona_old']['data']
    
    ax2 = plt.subplot(2,1,2)
    plt.plot_date(ftimes,np.sum(firepower,axis=(1,2)), '-k', label='Fire power waroona_run1')
    plt.plot_date(run1_xy[0], run1_xy[1], '--k', label='PFT waroona_run1')
    plt.plot_date(old_xy[0], old_xy[1], '--g', label='PFT waroona_old')
    
    plt.legend(loc='best')
    # ylabel and units
    plt.ylabel('Gigawatts')
    # format ticks
    #plt.xticks(hours[::4]) # just show one tick per 4 hours
    ax2.xaxis.set_major_locator(mdates.HourLocator(interval=4)) # every 4 hours
    ax2.xaxis.set_minor_locator(mdates.HourLocator())
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    fio.save_fig('mixed',_sn_,'firepower.png',plt)
### Run the stuff

#fireplan_waroona()
fire_power_waroona()