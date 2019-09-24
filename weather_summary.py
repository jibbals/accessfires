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
from datetime import datetime

from utilities import utils, plotting, fio

## GLOBAL
#Script name
_sn_ = 'weather_summary'

def plot_weather_summary(U,V,W, height, Q=None, 
                         model_version='waroona_oldold',
                         ext='.png',
                         datelimit=None,
                         dvi=150):
    '''
    Show horizontal slices of horizontal and vertical winds averaged between several vertical levels
    TODO: Show clouds summed over several levels if they are available
    
    INPUTS: 
        u: EW winds (cube [ t, levs, lats, lons ])
        v: NS winds
        w: Vert winds
        height: level heights for titles (array [ levs ])
    '''
    
    dtimes = utils.dates_from_iris(U)
    
    extentname=['sirivan','waroona'][dtimes[0].year==2016]
    row1 = (100<=height) * (height<500)
    row2 = (500<=height) * (height<1500)
    row3 = (1500<=height) * (height<3000)
    row4 = (3000<=height) * (height<5000)
    
    # vertical wind colourbar is constant
    wcmap=plotting._cmaps_['verticalvelocity']
    wnorm=colors.SymLogNorm(0.25) # linear to +- 0.25, then log scale
    wcontours=np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0]))
    
    for ti, dtime in enumerate(dtimes):

        if datelimit is not None:
            if dtime>datelimit:
                break
        
        
        fig = plt.figure(figsize=[10,10])
        for ii,row in enumerate([row1,row2,row3,row4]):
            plt.subplot(4,2,ii*2+1)
            
            Ui = U[ti,row]
            Vi = V[ti,row]
            
            # mean of levels # .data is masked array, just use base array
            Ur = np.mean(Ui.data.data,axis=0)
            Vr = np.mean(Vi.data.data,axis=0)
            windspeed = np.hypot(Ur,Vr)
            
            # value skip to vectors aren't so dense
            n_y,n_x = Ur.shape
            vsv = n_y // 16; vsu = n_x // 16 # about 16x16
            skip = (slice(None,None,vsv),slice(None,None,vsu))
            
            Y=Ui.coord('latitude').points
            X=Vi.coord('longitude').points
            Uvs,Vvs = Ur[skip], Vr[skip]
            
            # Normalize the arrows:
            
            Uvs = Uvs / np.sqrt(Uvs**2 + Vvs**2);
            Vvs = Vvs / np.sqrt(Uvs**2 + Vvs**2);
            
            plt.contourf(X, Y, windspeed, 10)
            
            plotting.map_add_locations_extent(extentname, hide_text=ii>0)
            
            plt.colorbar(ticklocation=ticker.MaxNLocator(5),pad=0)
            plt.quiver(X[::vsu], Y[::vsv], Uvs, Vvs, scale=30)
            plt.xticks([],[])
            plt.yticks([],[])
            plt.ylabel("%.0f - %.0f m"%(height[row][0], height[row][-1]))
            if ii==0:
                plt.title('horizontal winds (m/s)')
            
            # Now plot vertical motion
            plt.subplot(4,2,ii*2+2)
            Wi = W[ti,row]
            Wr = np.mean(Wi.data.data,axis=0)
            
            cs = plt.contourf(X,Y, Wr, wcontours, cmap=wcmap,norm=wnorm)
            plotting.map_add_locations_extent(extentname, hide_text=True)
            
            plt.xticks([],[])
            plt.yticks([],[])
            if ii==0:
                plt.title('vertical motion (m/s)')
        
        # add vert wind colourbar
        cbar_ax = fig.add_axes([0.905, 0.4, 0.01, 0.2])# X Y Width Height
        fig.colorbar(cs, cax=cbar_ax, format=ticker.ScalarFormatter(), pad=0)
        
        plt.suptitle("%s weather "%model_version + dtime.strftime("%Y %b %d %H:%M (UTC)"))
        dstamp = dtime.strftime("%Y%m%d%H%M")
        #pname = "figures/model_version/weather_summary/fig_%s%s"%(dstamp,ext)
        fio.save_fig(model_version,_sn_, dstamp, plt, ext=ext)

def read_and_plot_model_run(model_version='waroona_oldold',
                            dtimes=None):
    '''
    '''
    extentname = model_version.split('_')[0]
    extent = plotting._extents_[extentname]
    if dtimes is None:
        dtimes = fio.model_outputs['waroona_run1']['filedates']
    
    # font sizes etc
    plotting.init_plots()
    
    cubes = fio.read_model_run(model_version, fdtime=dtimes, extent=extent, 
                               add_winds=True)
    
    u,v,s = cubes.extract(['u','v','s'])
    w, = cubes.extract('upward_air_velocity')
    height = w.coord('level_height').points
    
    plot_weather_summary(u, v, w, height,
                         model_version=model_version)

if __name__=='__main__':
    dtimes = None
    if True: # for testing
        dtimes=[datetime(2016,1,6,8)]
    
    for mv in ['waroona_oldold','waroona_oldold','waroona_run1']:
        read_and_plot_model_run(mv, dtimes=dtimes)

