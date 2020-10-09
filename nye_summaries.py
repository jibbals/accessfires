#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue 20200609(YYYYMMDD)
    Weather summary figures - with edits so that NYE run plots properly
@author: jesse
"""
import matplotlib
matplotlib.use("Agg",warn=False)

from matplotlib import colors, ticker, patches
import matplotlib.dates as mdates

import numpy as np
import matplotlib.pyplot as plt
import iris # minimal fio done here 
import cartopy.crs as ccrs
import warnings
from datetime import datetime, timedelta

from utilities import utils, plotting, fio, constants
from weather_summary import plot_weather_summary, weather_series

## GLOBAL
#Script name
_sn_ = 'weather_summary'
__cloud_thresh__ = constants.cloud_threshold

def weather_summary_model(model_version='NYE_run1',
                          fdtimes=None,
                          zoom_in=None,
                          subdir=None,
                          HSkip=None,
                          hwind_limits=None):
    '''
    Read model run output hour by hour, running plot_weather_summary on each
    time slice. Can subset to just show fdtimes, and can also zoom to specific 
    lat lon box, which saves into the 'zoomed' subfolder
    '''
    # font sizes etc
    plotting.init_plots()
    
    extentname = model_version.split('_')[0]
    extent = constants.extents[extentname]
    if zoom_in is not None:
        extentname=None
        extent = zoom_in
    if fdtimes is None:
        fdtimes = fio.run_info[model_version]['filedates']
    FF = None
    
    # read one hour at a time, plot each available time slice
    for fdtime in fdtimes:
        
        # read cubes
        cubes = fio.read_model_run(model_version, fdtime=[fdtime], extent=extent, 
                                   add_winds=True,
                                   HSkip=HSkip)
        u,v = cubes.extract(['u','v'])
        w, = cubes.extract('upward_air_velocity')
        clouds = cubes.extract('qc')
        qc = None
        if len(clouds) == 1:
            qc = clouds[0]
        #print("DEBUG: lat coord", w.coord('latitude'))
        lat = w.coord('latitude').points
        lon = w.coord('longitude').points
        # convert lat/lon to distance...
        latm = np.linspace(0, utils.distance_between_points([lat[0],lon[0]],[lat[-1],lon[0]]), lat.size)
        lonm = np.linspace(0, utils.distance_between_points([lat[0],lon[0]],[lat[0],lon[-1]]), lon.size)

        ## LATS ARE NOT EXACTLY REGULAR?!!!
        print("DEBUG: pre-modify: lat0:",type(lat[0]), len(lat), lat[:3])
        #lat1 = np.around(lat,4) # This doesn't fix the lat array for some reason
        #print("DEBUG: post-around: lat1:",type(lat1[0]), len(lat1), lat1[:3], np.diff(lat1[:3]))
        lat2 = np.linspace(lat.min(), lat.max(), lat.size)
        print("DEBUG: post-recreate: lat2:",type(lat2[0]), len(lat2), lat2[:5], np.diff(lat2[:5]))
        height = w.coord('level_height').points
        dtimes = utils.dates_from_iris(u)
        # read fire front
        ff, = fio.read_fire(model_version, dtimes, extent=extent, HSkip=HSkip)
        
        # for each time slice create a weather summary plot
        for i,dtime in enumerate(dtimes):
            
            ui, vi = u[i].data.data, v[i].data.data
            wi = w[i].data.data
            qci = None
            if qc is not None:
                qci = qc[i].data.data
            
            plot_weather_summary(ui, vi, wi, height, lat2, lon, 
                                 extentname=extentname,
                                 Q = qci, FF=None, 
                                 hwind_limits=hwind_limits)
            
            plt.suptitle("%s weather "%model_version + dtime.strftime("%Y %b %d %H:%M (UTC)"))
            if zoom_in is not None and subdir is None: 
                subdir = 'zoomed'
            fio.save_fig(model_version, _sn_, dtime, plt, subdir=subdir)


if __name__=='__main__':
    mr = 'NYE_run1'
    NYE_hours = fio.run_info[mr]['filedates']
    NYE_zoom = [149.2,150.05, -37.5, -36.85]
    
    CanPlotTiff = True # set to True if not on NCI
    
    if False:
        ## NYE Weather summary figures:
        _sn_ = "weather_summary"
        weather_summary_model('NYE_run1',zoom_in=NYE_zoom,HSkip=None, 
                fdtimes=NYE_hours, hwind_limits=[0,30],)
    
    if CanPlotTiff:
        # Make a helper figure to show where the summary is located
        f,ax,proj = plotting.map_tiff_qgis("SEForests.tiff",)
        latlon_proj = ccrs.PlateCarree()
        ## Add box around zoomed in area
        xy = [NYE_zoom[0], NYE_zoom[2]]
        width = NYE_zoom[1]-NYE_zoom[0]
        height = NYE_zoom[3]-NYE_zoom[2]
        ax.add_patch(patches.Rectangle(xy=xy,
                                       width=width,
                                       height=height,
                                       #facecolor=None,
                                       fill=False,
                                       edgecolor='blue',
                                       linewidth=2,
                                       #linestyle='-',
                                       alpha=.7, 
                                       transform=latlon_proj))
        fio.save_fig("NYE_run1","weather_summary","weather_summary_extent.png",plt=plt)
    
    ## Now lets look at another thing!?

