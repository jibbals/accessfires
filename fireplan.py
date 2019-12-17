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
from matplotlib import patheffects
from matplotlib import ticker, colors, patches
import numpy as np
from datetime import datetime, timedelta
import iris

from cartopy import crs as ccrs

import PFT_work
from utilities import plotting, utils, fio

###
## GLOBALS
###
_sn_ = 'fireplan'
PFT = {'waroona_run1':{'data':np.array([56.5, 42.1, 34.6, 332.6, 53.2]), # manually calculated PFT
                       'units':'Gigawatts',
                       'time':np.array([datetime(2016,1,5,15,10),
                                        datetime(2016,1,6,5,10),
                                        datetime(2016,1,6,6,10),
                                        datetime(2016,1,6,8,10),
                                        datetime(2016,1,6,12,10)]), # calculated at these times in UTC
                       'latlon':plotting._latlons_['fire_waroona_upwind'], # ~ 1km from fire
                       'latlon_stamp':'fire_waroona_upwind',
                       'style':'--',
                       'color':'k',
                       },
       'waroona_old':{'data':np.array([61.4, 67.9, 176.4, 145.7]), # manually calculated PFT
                      'units':'Gigawatts',
                      'time':np.array([datetime(2016,1,5,15),
                                       datetime(2016,1,6,5,1),
                                       datetime(2016,1,6,6),
                                       datetime(2016,1,6,8,1),]), # calculated at these times in UTC
                      'latlon':plotting._latlons_['fire_waroona_upwind'], # ~ 1km from fire
                      'latlon_stamp':'fire_waroona_upwind',
                      'style':'--',
                      'color':'g',
                       },
       'sirivan_run1':{'data':None, # manually calculated PFT
                       'units':'Gigawatts',
                       'time':None, # calculated at these times in UTC
                       'latlon':plotting._latlons_['fire_sirivan_upwind'], # ~ 1km from fire
                       'latlon_stamp':'fire_sirivan_upwind',
                       'style':'--',
                       'color':'k',
                       },
      }

## Calculations for PFT:

## 201601060600OLD: q_ML=8.5, th_ML=34, z_fc=680mbar=3.5km, dth=6, U=8m/s :: PFT = .3*3.5**2*6*8 = 176.4 GW
## 201601060800OLD: q_ML=8, th_ML=35, z_fc=685mbar=3.4km, dth=4, U=10.5 :: PFT = .3*3.4**2*4*10.5 = 145.7 GW

## 201601060810RUN1: q_ML=8, th_ML=35.5, z_fc=690mbar=3.33km, dth=2, U=8 :: PFT = .3*3.33**2*2*8 = 53.2 GW

## 201601051500OLD: q_ML=8  , th_ML=35C, star_SP=(10  , 55C), z_fc=690mbar=3.24km, dth=1.5, U=13m/s ::: PFT = .3*3.24**2*1.5*13 = 61.4GW
## 201601060501OLD: q_ML=9  , th_ML=34C, star_SP=(11  , 54C), z_fc=710mbar=3.07km, dth=4.0, U= 6m/s ::: PFT = .3*3.07**2*4*6 = 67.9GW
## 201601051510NEW: q_ML=8  , th_ML=35C, star_SP=(10  , 55C), z_fc=710mbar=3.07km, dth=2.0, U=10m/s ::: PFT = .3*3.07**2*2*10 = 56.5GW

## 201601061210NEW: q_ML=7.2, th_ML=35.5, z_fc = 695mbar=4km, dth=4.5, U=15.4m/s ::: PFt = .3*4**2*4.5*15.4 = 332.6GW



def firepower_from_cube(shcube):
    """
    calculate and return firepower over time in GWatts
    
    Inputs
    ======
        shcube: sensible heat flux in Watts/m2
    """
    
    lon,lat = shcube.coord('longitude'), shcube.coord('latitude')
    ### get areas in m2
    # Add boundaries to grid
    if lat.bounds is None:
        lat.guess_bounds()
    if lon.bounds is None:
        lon.guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(shcube)

    firepower = shcube.data.data * grid_areas # W/m2 * m2
    return firepower/1e9 # Watts to Gigawatts

def pft_altitude_vs_pressure(model_run='waroona_run1', latlon=plotting._latlons_['fire_waroona_upwind'],
                             mbar_to_watch=700, datetimes=[datetime(2016,1,5,15)]):
    """
    Retrieve altitudes around 800-600 mbar over a specific lat lon over time.
    This is to calculate z_fc for the PFT calculation
    """
    ## First read the hourly z and
    extent = plotting._extents_[model_run.split('_')[0]]
    cubes = fio.read_model_run(model_run, fdtime=datetimes, extent=extent, add_topog=False,add_z=True)
    z,p = cubes.extract(['z_th','air_pressure'])
    p.convert_units('mbar')
    z.convert_units('km')
    cubetimes=utils.dates_from_iris(p)
    dstamp = cubetimes[0].strftime("%b %d %H:%M(UTC)")

    # pull out latlon that we are watching:
    lat,lon = latlon
    z0 = z.interpolate([('longitude',lon),('latitude',lat)],
                           iris.analysis.Linear())
    p0 = p.interpolate([('longitude',lon),('latitude',lat)],
                           iris.analysis.Linear())
    nt,nz = p0.shape

    z0,p0 = z0.data.data, p0.data.data

    pind = np.zeros(nt,dtype=np.int)
    # where is p0 closest to our watched pressure?
    for i in range(nt):
        pind[i] = np.argmin(np.abs(p0[i] - mbar_to_watch))
    # plot example scatter and z_fc grab

    plt.subplot(2,1,1)
    plt.scatter(p0[0], z0)
    plt.xlim([1000,500])
    plt.ylim([0,7])
    plt.xlabel('pressure [mbar]')
    plt.ylabel('altitude [km]')
    plt.plot([mbar_to_watch,mbar_to_watch],[0, z0[pind[0]]], color='red')
    plt.plot([1000, mbar_to_watch],[z0[pind[0]], z0[pind[0]]], color='red', label='closest to %d mbar'%mbar_to_watch)
    plt.legend(loc='best')
    plt.title('EG finding z$_{fc}$ at %s'%dstamp )


    plt.subplot(2,1,2)
    plt.plot_date(cubetimes,z0[pind])
    plt.title("z$_{fc}$",y=0.73)
    plt.ylabel('altitude [km]')

def fireplan(ff, fire_contour_map = 'autumn',
             show_cbar=True, cbar_XYWH= [0.2, 0.6, .2, .02],
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
        fig,... arguments to plotting.map_google()
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

    # First plot google map
    ## TODO: THIS PROCEDURE!!
    fig, gax, gproj = plotting.satellite_extent(extent=extent, fig=fig,
                                                subplot_row_col_n=subplot_row_col_n,
                                                draw_gridlines=False)
    
    # add waroona and waroona fire if we're over there
    latlons = plotting._latlons_
    if extent[0]<latlons['waroona'][1]: 
        plotting.map_add_nice_text(gax, 
                                   latlons=[latlons['waroona'],latlons['fire_waroona']], 
                                   texts=['Waroona',''],
                                   fontcolors='wheat',
                                   markers=['o','*'], 
                                   markercolors=['grey','red'])
    else:
        plotting.map_add_nice_text(gax,
                                   latlons=[latlons['uarbry'],latlons['cassillis'],latlons['fire_sirivan']],
                                   texts=['Uarbry','Cassillis',''],
                                   fontcolors='wheat',
                                   markers=['o','o','*'],
                                   markercolors=['k','k','red'],
                                   markersizes=4)
    
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

    return fig, gax, gproj

def fireplan_summary(model_run='waroona_run1', google=False):
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
             show_cbar=True, cbar_XYWH=[.2,.64,.2,.02],
             google=google)
    
    # Hourly for the rest of the stuff
    ftimes = utils.dates_from_iris(FFront)
    hourinds = [(ft.minute==0) and (ft.second==0) for ft in ftimes]
    FFront, SHeat, FSpeed = FFront[hourinds], SHeat[hourinds], FSpeed[hourinds]

    # plotting.map_add_locations(['fire_waroona_upwind'],text=['F160'],dx=-.001,dy=-.006, proj=ccrs.PlateCarree())
    ## attribute fire power from main part of fire
    #TODO

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

def fire_power_waroona():

    # Read fire output
    extentname = 'waroonaz' # fire zoomed
    #front_start=datetime(2016,1,6,1)
    extent = plotting._extents_[extentname]

    for model_run in ['waroona_run1','waroona_old']:
        # read all the fire data
        ff, sh = fio.read_fire(model_run=model_run, dtimes=None,
                               extent=extent, firefront=True, sensibleheat=True)
        lon,lat = ff.coord('longitude'), ff.coord('latitude')
        # just read hourly for fireplan
        ftimes = utils.dates_from_iris(ff)

        ## First plot fire front contours over google map in a 211 subplot
        fig = plt.figure(figsize=[8,8])

        _, ax, gproj = fireplan(ff,
                                fig=fig, subplot_row_col_n=[2,1,1],
                                show_cbar=True, cbar_XYWH=[.2,.64,.2,.02])
        
        plotting.map_add_nice_text(ax,
                                   latlons=[plotting._latlons_['fire_waroona_upwind']],
                                   texts='PFT',)
        #plotting.map_add_locations(['fire_waroona_upwind'], text=['PFT'],
        #                           dx=-.001, dy=-.006, proj=ccrs.PlateCarree())
        
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

        ax2 = plt.subplot(2,1,2)

        firepower = firepower_from_cube(sh)
        
        # Plot firepower
        plt.plot_date(ftimes,np.sum(firepower,axis=(1,2)), '-'+PFT[model_run]['color'], label='Fire power '+model_run)

        # also PFT values:
        for mr in ['waroona_run1','waroona_old']:
            linestyle = PFT[mr]['style']+PFT[mr]['color']
            run_x,run_y = PFT[mr]['time'],PFT[mr]['data']
            plt.plot_date(run_x, run_y, linestyle, label='PFT '+mr)

        plt.legend(loc='best')
        # ylabel and units
        plt.ylabel('Gigawatts')
        # format ticks
        #plt.xticks(hours[::4]) # just show one tick per 4 hours
        ax2.xaxis.set_major_locator(mdates.HourLocator(interval=4)) # every 4 hours
        ax2.xaxis.set_minor_locator(mdates.HourLocator())
        ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        fio.save_fig(model_run, _sn_, 'firepower.png',plt)

def firepower_comparison(runs=['waroona_new1','waroona_old'], localtime=False):
    """
    Plot overlaid time series of two model runs fire power from integral of intensity
    """
    lat,lon = PFT[runs[0]]['latlon']
    fpextent = plotting._extents_[runs[0].split('_')[0]+'z']
    pftextent = [lon-.01, lon+.01, lat-.01, lat+.01]
    offset = timedelta(hours=8)
    if lon > 130:
        offset = timedelta(hours=10)
    plt.figure(figsize=[10,5]) # time series
    
    labels = runs
    #labels = ['new','orig'] # just for poster
    
    for i,run in enumerate(runs):
        # read all the fire data
        sh, = fio.read_fire(model_run=run, dtimes=None,
                            extent=fpextent, firefront=False, sensibleheat=True)
        
        ftimes = utils.dates_from_iris(sh)
        if localtime:
            ftimes = np.array([ft + offset for ft in ftimes ])
        firepower = np.sum(firepower_from_cube(sh), axis=(1,2)) # over time
        
        # remove zeros:
        prefire = np.isclose(np.cumsum(firepower), 0)
        firepower[prefire] = np.NaN
        
        # Plot firepower
        plt.plot_date(ftimes, firepower, '-'+PFT[run]['color'], \
                      label=labels[i])
        
        # also PFT values:
        linestyle = PFT[run]['style']+PFT[run]['color']
        # calculate using kevin's code
        print("DEBUG: reading sirivan")
        cubes = fio.read_model_run(run, extent=pftextent,
                                   add_theta=True, add_winds=True, add_RH=True, 
                                   add_z=True, add_topog=True)
        
        print("DEBUG:",cubes)
        utimes = utils.dates_from_iris(cubes.extract('u')[0])
        if localtime:
            utimes = np.array([ut + offset for ut in utimes ])
        
        pft = PFT_work.PFT_from_cubelist(cubes,latlon=[lat,lon])
        plt.plot_date(utimes, pft, linestyle, label='PFT '+labels[i])
    
    # fix up labels, axes
    plt.legend(loc='best')
    # ylabel and units
    plt.ylabel('Gigawatts')
    plt.title('Firepower')
    # format ticks
    #plt.xticks(hours[::4]) # just show one tick per 4 hours
    ax = plt.gca()
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=4)) # every 4 hours
    ax.xaxis.set_minor_locator(mdates.HourLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    fio.save_fig(run, _sn_, 'firepower_comparison.png',plt)
    

if __name__=='__main__':
    ### Run the stuff
    
    
    ## create firepower time series
    #firepower_comparison(runs=['waroona_run1'])
    #firepower_comparison(runs=['waroona_old'])
    #firepower_comparison(runs=['waroona_run1','waroona_old'])
    
    #firepower_comparison(runs=['sirivan_run1'])
    
    for mr in ['sirivan_run1','waroona_run1','waroona_old','waroona_run2']:
        # zoomed extent
        extent = plotting._extents_[mr.split('_')[0]+'z'] 
        # read all the fire data
        ff, = fio.read_fire(model_run=mr, dtimes=None,
                            extent=extent, firefront=True)
        
        # first plot just the fireplan on it's own
        fireplan(ff, show_cbar=True, cbar_XYWH=[.2,.24,.2,.02])
        fio.save_fig(mr,_sn_,'fire_outline.png',plt,dpi=300)
        
        # Now run summary figures
        fireplan_summary(model_run=mr)
    
    