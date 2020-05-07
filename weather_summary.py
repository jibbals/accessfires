#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 08:45:14 2019
    Weather summary looking at winds and clouds
@author: jesse
"""
import matplotlib
matplotlib.use("Agg",warn=False)

from matplotlib import colors, ticker

import numpy as np
import matplotlib.pyplot as plt
import warnings
from datetime import datetime

from utilities import utils, plotting, fio, constants

## GLOBAL
#Script name
_sn_ = 'weather_summary'
__cloud_thresh__ = constants.cloud_threshold

def plot_weather_summary(U,V,W, height, lat, lon, extentname, 
                         Q=None, FF=None, Streamplot=True):
    '''
    Show horizontal slices of horizontal and vertical winds averaged between 
    several vertical levels. Also shows clouds (Q) and fires (FF) with contour outline.
    
    INPUTS: 
        U,V,W: wind speed in lon, lat, vert dimension m/s [z,lat,lon]
        height: altitude ASL or AGL m [z,lat,lon] or [z]
        height: level heights for titles (array [ levs ])
    '''
    
    row1 = (100<=height) * (height<500)
    row2 = (500<=height) * (height<1500)
    row3 = (1500<=height) * (height<3000)
    row4 = (3000<=height) * (height<5000)
    # todo row 5
    row5 = (5000<height) * (height<10000)
    
    # vertical wind colourbar is constant
    wcmap=plotting._cmaps_['verticalvelocity']
    wnorm=colors.SymLogNorm(0.25) # linear to +- 0.25, then log scale
    wcontours=np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0]))
    
    fig = plt.figure(figsize=[10,10])
    for ii,row in enumerate([row1, row2, row3, row4, row5]):
        plt.subplot(5,2,ii*2+1)
        
        Ui = U[row]
        Vi = V[row]
        Si = np.hypot(Ui,Vi)
        
        # mean of levels # .data is masked array, just use base array
        Ur = np.mean(Ui, axis=0)
        Vr = np.mean(Vi, axis=0)
        Sr = np.mean(Si, axis=0)
        
        # wind speed contourf
        plt.contourf(lon, lat, Sr, 10)
        plt.colorbar(ticklocation=ticker.MaxNLocator(5),pad=0)
        
        if extentname is not None:
            plotting.map_add_locations_extent(extentname, hide_text=True)
        
        if Streamplot:
            # how dense is streamplot
            density_x = 0.6
            density_y = 0.5
            #streamplot for horizontal winds
            plt.streamplot(lon,lat,Ur,Vr, 
                           color='k', 
                           density=(density_x, density_y))
        else:
            # value skip so vectors aren't so dense
            n_y, n_x = Ur.shape
            vsv = n_y // 16; vsu = n_x // 16 # about 16x16
            skip = (slice(None,None,vsv),slice(None,None,vsu))    
            Uvs,Vvs = Ur[skip], Vr[skip]
        
            # Normalize the arrows:
            Uvs = Uvs / np.sqrt(Uvs**2 + Vvs**2);
            Vvs = Vvs / np.sqrt(Uvs**2 + Vvs**2);
            
            # quiver plot
            plt.quiver(lon[::vsu], lat[::vsv], Uvs, Vvs, scale=30)
        
        # remove x and y ticks
        plt.xticks([],[])
        plt.yticks([],[])
        # add ylabel
        plt.ylabel("%.0f - %.0f m"%(height[row][0], height[row][-1]),fontsize=13)
        if ii==0:
            plt.title('horizontal winds (m/s)')
        
        # Now plot vertical motion
        plt.subplot(5,2,ii*2+2)
        Wi = W[row]
        Wr = np.mean(Wi,axis=0)
        
        cs = plt.contourf(lon, lat, Wr, wcontours, cmap=wcmap, norm=wnorm)
        if extentname is not None:
            plotting.map_add_locations_extent(extentname, hide_text=ii>0)
        
        # add cloud hatching
        if Q is not None:
            Qmax = np.max(Q[row], axis=0)
            # add hatches over where Qmax is greater than cloud threshhold
            plt.contourf(lon,lat,Qmax, [0,__cloud_thresh__,100], 
                         colors=['None']*3, hatches=[None,'/','//'],
                         extend='both',)
        
        if FF is not None:
            plotting.map_fire(FF, lat, lon)
        
        plt.xticks([],[])
        plt.yticks([],[])
        if ii==0:
            plt.title('vertical motion (m/s)')
        
    # reduce vert gap between subplots
    fig.subplots_adjust(hspace=0.1)
    # add vert wind colourbar
    cbar_ax = fig.add_axes([0.905, 0.4, 0.01, 0.2])# X Y Width Height
    fig.colorbar(cs, cax=cbar_ax, format=ticker.ScalarFormatter(), pad=0)
        

def weather_summary_model(model_version='waroona_run1',
                          fdtimes=None,
                          zoom_in=None,
                          HSkip=None):
    '''
    Read model run output hour by hour, running plot_weather_summary on each
    time slice. Can subset to just show fdtimes, and can also zoom to specific 
    lat lon box, which saves into the 'zoomed' subfolder
    '''
    # font sizes etc
    plotting.init_plots()
    
    extentname = model_version.split('_')[0]
    extent = plotting._extents_[extentname]
    if zoom_in is not None:
        extentname=None
        extent = zoom_in
    if fdtimes is None:
        fdtimes = fio.model_outputs[model_version]['filedates']
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
        lat = w.coord('latitude').points
        lon = w.coord('longitude').points
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
            if ff is not None:
                FF = ff[i].data.data
            
            plot_weather_summary(ui, vi, wi, height, lat, lon, 
                                 extentname=extentname,
                                 Q = qci, FF=FF)
            
            plt.suptitle("%s weather "%model_version + dtime.strftime("%Y %b %d %H:%M (UTC)"))
            subdir=None
            if zoom_in is not None: 
                subdir = 'zoomed'
            fio.save_fig(model_version,_sn_, dtime, plt, subdir=subdir)


if __name__=='__main__':
    
    ## run for all of waroona_run2 datetimes
    weather_summary_model(model_version='waroona_run3')
    
    ## run zoomed in
    #zoom_in = plotting._extents_['sirivans']
    #zoom_in = None # or not
    #weather_summary_model('waroona_run3',zoom_in=zoom_in,HSkip=None, 
    #                      fdtimes=[datetime(2016,1,6,13)])

    print("INFO: weather_summary.py done")

