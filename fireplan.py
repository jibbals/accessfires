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
             show_cbar=True, cbar_XYWH= [0.2,0.6,.3,.02],
             fig=None,subplot_row_col_n=None,draw_gridlines=False,gridlines=None):
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
    extent = [np.min(lon),np.max(lon), np.min(lat),np.max(lat)] 
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
                                           fig=fig, subplot_row_col_n=subplot_row_col_n, 
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

def fireplan_summary(model_run='waroona_run1'):
    '''
    Show fire outline over time at waroona
    '''
    # Read fire output
    extentname1 = model_run.split('_')[0]
    extentname = extentname1+'z' # fire zoomed
    
    extent = plotting._extents_[extentname] 
    
    # just read hourly
    hourly = [datetime(2016,1,5,15) + timedelta(hours=x) for x in range(24)]
    FFront, SHeat, FSpeed = fio.read_fire(model_run=model_run ,dtimes=hourly, extent=extent,
                                          firefront=True, sensibleheat=True, firespeed=True)
    lon,lat = FFront.coord('longitude').points, FFront.coord('latitude').points
        
    ## PLOTTING
    
    ## First plot fire front contours over google map in a 211 subplot
    fig = plt.figure(figsize=[8,9])
    
    fireplan(FFront, extentname=extentname,
             fig=fig, subplot_row_col_n=[2,1,1],
             show_cbar=True, cbar_XYWH=[.2,.6,.3,.02])
    
    # plotting.map_add_locations(['fire_waroona_upwind'],text=['F160'],dx=-.001,dy=-.006, proj=ccrs.PlateCarree())
    
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
    cs=plt.contourf(lon, lat, maxspeed.T, 30, 
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
        hourinds = [ft.minute==0 for ft in ftimes]
        hours = ftimes[hourinds]
        
        ## First plot fire front contours over google map in a 211 subplot
        fig = plt.figure(figsize=[8,8])
        
        _, ax, gproj = fireplan(ff[hourinds],
                                fig=fig, subplot_row_col_n=[2,1,1],
                                show_cbar=True, cbar_XYWH=[.2,.6,.3,.02])
        
        plotting.map_add_locations(['fire_waroona_upwind'],text=['PFT'],dx=-.001,dy=-.006, proj=ccrs.PlateCarree())
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
        
        
        ### get areas in m2
        # Add boundaries to grid
        lat.guess_bounds()
        lon.guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(ff)
        
        firepower = sh.data.data * grid_areas # W/m2 * m2
        firepower = firepower/1e9 # Watts to Gigawatts
        
        # Plot firepower
        plt.plot_date(ftimes,np.sum(firepower,axis=(1,2)), '-k', label='Fire power waroona_run1')
        
        # also PFT values:
        for mr, linestyle in zip(['waroona_run1','waroona_old'],['--k','--g']):
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

if __name__=='__main__':
    ### Run the stuff
    
    fireplan_summary()
    fire_power_waroona()