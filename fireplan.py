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
from matplotlib import patheffects, colors, image, ticker
import numpy as np
# for legend creation:
from matplotlib.lines import Line2D
from datetime import datetime,timedelta



import iris

from cartopy import crs as ccrs

from utilities import plotting, utils, fio

###
## GLOBALS
###
_sn_ = 'fireplan'


def fireplan(ff, fire_contour_map = 'autumn',
             show_cbar=True, cbar_XYWH= [0.65, 0.63, .2, .02],
             **kwtiffargs):
#             fig=None,subplot_row_col_n=None,
#             draw_grid=False,gridlines=None):
    '''
    show satellite map of extent, with fire front over time overplotted

    ARGUMENTS:
        ff: iris.cube.Cube with time, longitude, latitude dimensions
        fire_contour_map: how will firefront contours be coloured
        show_cbar: bool
            draw a little colour bar showing date range of contours?
        cbar_XYHW: 4-length list where to put cbar
        fig,... arguments to plotting.map_tiff_qgis()
    '''
    lon,lat = ff.coord('longitude').points, ff.coord('latitude').points
    _,ny,nx = ff.shape
    extent = [np.min(lon),np.max(lon), np.min(lat),np.max(lat)]

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
    extentname='sirivan' if extent[0]>130 else 'waroonaf'
    fig, ax = plotting.map_tiff_qgis(
        fname=extentname+'.tiff',
        extent=extent, 
        **kwtiffargs,
        )
    if 'waroona' in extentname:
        plotting.map_add_nice_text(ax,[plotting._latlons_['waroona'],plotting._latlons_['yarloop']],
                                   texts=['Waroona','Yarloop'], fontsizes=14)
    else:
        plotting.map_add_nice_text(ax,[plotting._latlons_['uarbry']],
                                   texts=['Uarbry'], fontsizes=14)
    
    # plot contours at each hour
    tstamp=[]
    for ii,dt in enumerate(ftimes[hourinds]):
        offset = 8 if "waroona" in extentname else 10
        LT = dt + timedelta(hours=offset)
        color=[rgba[ii]]
        # Colour waroona by day of firefront
        if "waroona" in extentname:
            dnum = int(LT.strftime("%d"))
            if dnum == 6:
                color = 'red'
            elif dnum==7:
                color="orange"
            else:
                color = "grey"
                
        tstamp.append(LT.strftime('%b %d, %H%m(local)'))

        ffdata = ff_f[hourinds][ii].data.data
        
        linewidth=1 + (ii == nt-1) # make final hour thicker
        fire_line = plt.contour(lon, lat, ffdata, np.array([0]),
                                colors=color, linewidths=linewidth,
                                )
    
    # final contour gets bonus blue line
    final_line=plt.contour(lon,lat,ff_f[-1].data.data, np.array([0]), 
                           linestyles='dotted',
                           colors='cyan', linewidths=1)
    clbls = plt.clabel(final_line,[0],fmt=LT.strftime('%H:%M'), 
                       inline=True, colors='wheat')
    plt.setp(clbls, path_effects=[patheffects.withStroke(linewidth=3, foreground="k")])

    ## Add tiny colour bar showing overall time of fire
    if show_cbar:
        # create an axis somewhere to add a colorbar
        cax = fig.add_axes(cbar_XYWH)
        #cax.axes.get_xaxis().set_visible(False)
        cax.axes.get_yaxis().set_visible(False) # horizontal alignment, no yax
        cbar = matplotlib.colorbar.ColorbarBase(cax,cmap=plt.get_cmap(fire_contour_map), orientation='horizontal')
        cbar.set_ticks([0,1])
        cbar.set_ticklabels([tstamp[0],tstamp[-1]])
        # make xtick labels have black outline on wheat text color
        cbxtick_obj = plt.getp(cax, 'xticklabels') # get xtick labels
        plt.setp(cbxtick_obj, color='wheat',
                 path_effects=[patheffects.withStroke(linewidth=3, foreground="k")])
        # return focus to newly created plot
        plt.sca(ax)
    plotting.scale_bar(ax, ccrs.PlateCarree(), 10)

    return fig, ax

def fireplan_summary(model_run='waroona_run1',
                     day1=True,day2=True,
                     just_fireplan=False):
    '''
    Show fire outline over time at waroona
    '''
    # Read fire output
    extentname1 = model_run.split('_')[0]
    extentname = extentname1+'f' # area affected by fire
    extent = plotting._extents_[extentname]
    
    if 'path_firefront2' not in fio.model_outputs[model_run].keys():
        day2 = False
    FFront, SHeat, FSpeed = fio.read_fire(model_run=model_run ,dtimes=None, 
                                          extent=extent, firefront=True, 
                                          sensibleheat=True, firespeed=True,
                                          day1=day1, day2=day2)
    lon,lat = FFront.coord('longitude').points, FFront.coord('latitude').points

    ## PLOTTING: first just show firefronts hourly contour
    fig,ax = fireplan(FFront, show_cbar=True, cbar_XYWH=[0.18,0.075,.2,.02])
    fio.save_fig(model_run,_sn_,"fireplan",plt)

    if just_fireplan:
        return None
    
    ## Bigger figure starts with plot fire front contours
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
    
    fio.save_fig(model_run, _sn_, 'firesummary', plt)
    
    
def fireplan_comparison(model_runs=['waroona_old','waroona_run1','waroona_run2','waroona_run3'],
                        colors = ['red','orange','teal','magenta'],
                        extent=plotting._extents_['waroona'],
                        mapname='waroona.tiff',
                        figname='waroona_firespread'):
    """
    compare fire spread between multiple runs
    """
    
    # first pull tiff image into memory
    fig = plt.figure(figsize=[14,10])
    fig,ax = plotting.map_tiff_qgis(fname=mapname,extent=extent,fig=fig)
    # firefront coord system is lats and lons
    legend = []
    
    # for each model run, extract last fire front
    for mr, color in zip(model_runs, colors):
        FF, = fio.read_fire(model_run=mr,
                            dtimes=None,
                            extent=extent,
                            firefront=True,)
        lon,lat=FF.coord('longitude').points,FF.coord('latitude').points
        print("DEBUG:",FF.summary(shorten=True))
        # last firefront
        FFlast = FF[-1].data
        # first:
        firstind = np.where(np.min(FF.data,axis=(1,2))<0)[0][0]
        FFfirst = FF[firstind].data
        # contour the firefronts
        plt.contour(lon, lat, FFfirst, np.array([0]),
                    colors=color, 
                    linewidths=2,
                    alpha=0.5,
                    )
        
        plt.contour(lon, lat, FFlast, np.array([0]),
                    colors=color, 
                    linewidths=2,
                    alpha=0.75,
                    )
        legend.append(Line2D([0],[0],color=color,lw=2))
    
    if 'waroona' in model_runs[0]:
        plotting.map_add_nice_text(ax,
                                   [plotting._latlons_['waroona'],plotting._latlons_['yarloop']],
                                   texts=['Waroona','Yarloop'], fontsizes=14)
    
    # Make/Add legend
    ax.legend(legend, model_runs)
    plt.tight_layout()
    fio.save_fig('project', _sn_, figname, plt)

def fireplan_vs_isochrones():
    """
    Show isochrones over waroona, look at flux over yarloop and compare to iso
    4 panels: 
        A: isochrones, 
        B: FF contours (hourly coloured by localtime day)
        C: flux at 3pm on day2? zoomed into yarloop a bit
        D: flux at 8pm on day2? zoomed a bit
    """
    mr='waroona_run3'
    fig = plt.figure(figsize=[14,12])
    ax1 = plt.subplot(2,2,1)
    ## Top left is the isochrones picture
    iso = image.imread('data/Waroona_Fire_Isochrones.png')
    ax1.imshow(iso)
    plt.xticks([],[])
    plt.yticks([],[])

    # Read fire output
    # area affected by fire
    extentA = [115.6,116.21, -33.05,-32.8]
    # area to zoom in on for yarloop
    extentB = [115.7,116.05, -33,-32.82]
    
    FFront, SHeat, = fio.read_fire(model_run=mr, dtimes=None, 
                                   extent=extentA, 
                                   firefront=True, 
                                   sensibleheat=True,
                                   day1=True, day2=True
                                   )
    lon, lat = FFront.coord('longitude').points, FFront.coord('latitude').points
    fdates=utils.dates_from_iris(FFront)

    ## PANEL B: firefronts hourly contours
    _,ax2 = fireplan(FFront, show_cbar=False, fig=fig, subplot_row_col_n=[2,2,2])
    
    ## PANEL C: flux at 3PM local time over yarloop (day 2?)
    dt1 = datetime(2016,1,7,7,30)
    dt2 = datetime(2016,1,7,12,30)
    for i,dt in enumerate([dt1,dt2]):
        dti = utils.nearest_date_index(dt,fdates)
        lt = dt + timedelta(hours=8)

        _,ax3,_ = plotting.map_tiff_qgis(
            "waroonaf.tiff", 
            extent=extentB, 
            fig=fig, 
            subplot_row_col_n=[2,2,3+i],
            locnames=['yarloop','waroona'],
            )
    
        cs,cbar = plotting.map_sensibleheat(SHeat[dti].data.data,
                                            lat,lon,
                                            colorbar=True,
                                            alpha=0.8)
        plt.title(
            lt.strftime('Sensible heat flux, %d %H:%M (LT)'),
            fontsize=14,
            #,y=.74)
            )
        plt.xticks([],[])
        plt.yticks([],[])
    plt.subplots_adjust(left=.03, right=.97, top=.96, bottom=.04,
            wspace=.03, hspace=.01)
    
    fio.save_fig(mr,_sn_,"fireplan_vs_isochrones.png",plt)


if __name__=='__main__':
    ### Run the stuff
    
    ##fireplan comparison
    if False:
        fireplan_comparison()
    
    if True:
        fireplan_vs_isochrones()
    
    ## Just create a fireplan figure:
    if False:
        fireplanruns = ['sirivan_run2_hr','waroona_old','waroona_run1','sirivan_run1','waroona_run3','waroona_run2']
        
        si_r2_hr = 'sirivan_run2_hr'
        si_r1 = 'sirivan_run1'
        
        extent = plotting._extents_['sirivans'] # synoptic extent
    
        ## Plot fireplan for high res run
        # read all the fire data
        ff, = fio.read_fire(model_run=si_r2_hr, dtimes=None,
                            extent=extent, firefront=True,
                            HSkip=5)
        
        fig,ax = fireplan(ff, show_cbar=True, cbar_XYWH=[.18,.3,.2,.02])
        fio.save_fig(si_r2_hr, _sn_, 'fireplan_hr.png', plt)
        
        ## Plot fireplan for sirivan original run
        # read all the fire data
        ff, = fio.read_fire(model_run=si_r1, dtimes=None,
                            extent=extent, firefront=True,
                            HSkip=5)
        
        # first plot just the fireplan on it's own
        fig,ax = fireplan(ff, show_cbar=True, cbar_XYWH=[.18,.3,.2,.02])
        fio.save_fig(si_r1, _sn_, 'fireplan_hr.png', plt)

    if False:
        mr = "waroona_run3"
        fireplan_summary(model_run=mr, day2=True, just_fireplan=True)
    
    
