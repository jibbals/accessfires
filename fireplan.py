#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 16:28:18 2019
    Show fire spread and intensity over time
@author: jesse
"""


import matplotlib
matplotlib.use('Agg',warn=False)
import matplotlib.pyplot as plt
from matplotlib import patheffects, colors
import numpy as np

import iris

from cartopy import crs as ccrs

from utilities import plotting, utils, fio

###
## GLOBALS
###
_sn_ = 'fireplan'


def fireplan(ff, fire_contour_map = 'autumn',
             show_cbar=True, cbar_XYWH= [0.65, 0.63, .2, .02],
             fig=None,subplot_row_col_n=None,
             draw_gridlines=False,gridlines=None):
    '''
    show satellite map of extent, with fire front over time overplotted

    ARGUMENTS:
        ff: iris.cube.Cube with time, longitude, latitude dimensions
        extentname: 'waroona' or 'sirivan'
        fire_contour_map: how will firefront contours be coloured
        show_cbar: bool
            draw a little colour bar showing date range of contours?
        cbar_XYHW: 4-length list where to put cbar
        fig,... arguments to plotting.map_tiff()
    '''
    lon,lat = ff.coord('longitude').points, ff.coord('latitude').points
    _,nx,ny = ff.shape
    extent = [np.min(lon),np.max(lon), np.min(lat),np.max(lat)]
    crs_data = ccrs.PlateCarree()

    # Get datetimes from firefront cube
    ftimes = utils.dates_from_iris(ff)
    
    # How many fire contours do we have?
    minff = np.min(ff.data,axis=(1,2))
    subset_with_fire = minff < 0
    ff_f = ff[subset_with_fire]
    ftimes = ftimes[subset_with_fire]
    
    # just read hourly
    hourinds = [(ft.minute==0) and (ft.second==0) for ft in ftimes]
    nt = np.sum(hourinds)
    
    ## fire contour colour map
    cmap = matplotlib.cm.get_cmap(fire_contour_map)
    rgba = cmap(np.linspace(0,1,nt))

    ## PLOTTING

    # First show satellite image and locations
    locname='waroona_big'
    if extent[0]>130: locname='sirivan'
    fig, gax, gproj = plotting.map_tiff(
        locname=locname,
        extent=extent, 
        fig=fig,
        subplot_row_col_n=subplot_row_col_n,
        show_grid=False
    )
    if 'waroona' in locname:
        plotting.map_add_nice_text(gax,[plotting._latlons_['waroona'],plotting._latlons_['yarloop']],
                                   texts=['Waroona','Yarloop'], fontsizes=14)
    else:
        plotting.map_add_nice_text(gax,[plotting._latlons_['uarbry']],
                                   texts=['Uarbry'], fontsizes=14)
    
    # plot contours at each hour
    utcstamp=[]
    for ii,dt in enumerate(ftimes[hourinds]):

        utcstamp.append(dt.strftime('%b %d, %HH(UTC)'))

        ffdata = ff_f[hourinds][ii].data.data
        
        linewidth=1 + (ii == len(hourinds)-1) # make final hour thicker
        plt.contour(lon, lat, ffdata.T, np.array([0]),
                    colors=[rgba[ii]], linewidths=linewidth,
                    transform=crs_data)
        
    # final outline contour
    final_line=plt.contour(lon,lat,ff_f[-1].data.data.T, np.array([0]), 
                           linestyles='dotted',
                           colors='cyan', linewidths=1, transform=crs_data)
    clbls = plt.clabel(final_line,[0],fmt=ftimes[-1].strftime('%H:%M'), inline=True, color='wheat')
    plt.setp(clbls, path_effects=[patheffects.withStroke(linewidth=3, foreground="k")])

    ## Add tiny colour bar showing overall time of fire
    if show_cbar:
        # create an axis somewhere to add a colorbar
        cax = fig.add_axes(cbar_XYWH)
        #cax.axes.get_xaxis().set_visible(False)
        cax.axes.get_yaxis().set_visible(False) # horizontal alignment, no yax
        cbar = matplotlib.colorbar.ColorbarBase(cax,cmap=plt.get_cmap(fire_contour_map), orientation='horizontal')
        cbar.set_ticks([0,1])
        cbar.set_ticklabels([utcstamp[0],utcstamp[-1]])
        # make xtick labels have black outline on wheat text color
        cbxtick_obj = plt.getp(cax, 'xticklabels') # get xtick labels
        plt.setp(cbxtick_obj, color='wheat',
                 path_effects=[patheffects.withStroke(linewidth=3, foreground="k")])
        # return focus to newly created plot
        plt.sca(gax)
    plotting.scale_bar(gax, gproj, 10)

    return fig, gax, gproj

def fireplan_summary(model_run='waroona_run1'):
    '''
    Show fire outline over time at waroona
    ARGUMENTS:
        main_box=[EWSN] extent of main fire, optionally attribute fire power from a subset of extent
    '''
    # Read fire output
    extentname1 = model_run.split('_')[0]
    extentname = extentname1+'z' # fire zoomed

    extent = plotting._extents_[extentname]

    FFront, SHeat, FSpeed = fio.read_fire(model_run=model_run ,dtimes=None, extent=extent,
                                          firefront=True, sensibleheat=True, firespeed=True)
    lon,lat = FFront.coord('longitude').points, FFront.coord('latitude').points

    ## PLOTTING

    ## First plot fire front contours over google map in a 211 subplot
    fig = plt.figure(figsize=[8,9])

    fireplan(FFront,
             fig=fig, subplot_row_col_n=[2,1,1],
             show_cbar=True, cbar_XYWH=[.2,.64,.2,.02])
    
    # Hourly for the rest of the stuff
    ftimes = utils.dates_from_iris(FFront)
    hourinds = [(ft.minute==0) and (ft.second==0) for ft in ftimes]
    FFront, SHeat, FSpeed = FFront[hourinds], SHeat[hourinds], FSpeed[hourinds]

    ## subplot 2
    ax2 = plt.subplot(4,1,3)
    maxflux = np.sum(SHeat.data.data,axis=0) + 0.01 # get rid of zeros
    levels = np.sort(np.union1d(np.power(10,np.arange(2,6)),5*np.power(10,np.arange(2,6))))
    cs = plt.contourf(lon, lat, maxflux.T,
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
    cs = plt.contourf(lon, lat, maxspeed.T, 30,
                    cmap=plotting._cmaps_['windspeed'])
    plt.colorbar(cs, orientation='horizontal', pad=0)
    plt.title('Max Firespeed (m/s ?)', y=.74)
    plt.xticks([],[])
    plt.yticks([],[])
    
    fio.save_fig(model_run, _sn_, 'fire_spread', plt)
    

if __name__=='__main__':
    ### Run the stuff
    
    ## Just create a fireplan figure:
    for mr in ['waroona_run3','waroona_run2']:
        #extent = plotting._extents_[mr.split('_')[0]+'z']  # zoomed waroona extent
        extent = [115.6, 116.2, -33.06, -32.8]
        # read all the fire data
        ff, = fio.read_fire(model_run=mr, dtimes=None,
                            extent=extent, firefront=True)
    
        # first plot just the fireplan on it's own
        fig,ax,proj = fireplan(ff, show_cbar=True, cbar_XYWH=[.18,.3,.2,.02])
        fio.save_fig(mr, _sn_, 'fireplan.png', plt)
    
    ## run fireplan and summary for all runs
    #for mr in ['waroona_run2','sirivan_run1','waroona_run1','waroona_old']:
    #    # zoomed extent
    #    extent = plotting._extents_[mr.split('_')[0]+'z'] 
    #    # read all the fire data
    #    ff, = fio.read_fire(model_run=mr, dtimes=None,
    #                        extent=extent, firefront=True)
    #    
    #    # first plot just the fireplan on it's own
    #    fireplan(ff, show_cbar=True, cbar_XYWH=[.2,.24,.2,.02])
    #    fio.save_fig(mr,_sn_,'fire_outline.png',plt,dpi=300)
    #    
    #    # Now run summary figures
    #    fireplan_summary(model_run=mr)
    
    
