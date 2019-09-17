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

import iris
import numpy as np
import matplotlib.pyplot as plt

from utilities import utils, plotting, fio


def plot_weather_summary(U,V,W, height, Q=None, 
                         model_version='waroona_oldold',
                         ext='.png', timedim_name='time',
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
    
    cubetimes = U.coord(timedim_name)
    dtimes = utils.dates_from_iris(cubetimes)
    
    extentname=['sirivan','waroona'][dtimes[0].year==2016]
    row1 = (100<=height) * (height<500)
    row2 = (500<=height) * (height<1500)
    row3 = (1500<=height) * (height<3000)
    row4 = (3000<=height) * (height<5000)
    
    for ti, dtime in enumerate(dtimes):

        if datelimit is not None:
            if dtime>datelimit:
                break
        
        plt.figure(figsize=[10,10])
        for ii,row in enumerate([row1,row2,row3,row4]):
            plt.subplot(4,2,ii*2+1)
            
            Ui = U[ti,row]
            Vi = V[ti,row]
            
            # mean of levels
            Ur = np.mean(Ui.data,axis=0)
            Vr = np.mean(Vi.data,axis=0)
            windspeed = np.hypot(Ur,Vr)
            
            # value skip to vectors aren't so dense
            n_y,n_x = Ur.shape
            vsv = n_y // 16; vsu = n_x // 16 # about 16x16
            skip = (slice(None,None,vsv),slice(None,None,vsu))
            
            Y=Ui.coord('latitude').points
            X=Vi.coord('longitude').points
            Uvs,Vvs = Ur[skip].data, Vr[skip].data
            
            # Normalize the arrows:
            
            Uvs = Uvs / np.sqrt(Uvs**2 + Vvs**2);
            Vvs = Vvs / np.sqrt(Uvs**2 + Vvs**2);
            
            plt.contourf(X, Y, windspeed.data)
            
            plotting.add_map_locations(extentname, hide_text=ii>0)
            
            plt.colorbar(pad=0)
            plt.quiver(X[::vsu], Y[::vsv], Uvs, Vvs, scale=30)
            plt.xticks([],[])
            plt.yticks([],[])
            plt.ylabel("%.0f - %.0f m"%(height[row][0], height[row][-1]))
            if ii==0:
                plt.title('horizontal winds (m/s)')
            
            # Now plot vertical motion
            plt.subplot(4,2,ii*2+2)
            Wi = W[ti,row]
            Wr = np.mean(Wi.data,axis=0)
            
            # colour bar stuff
            cmap=plotting._cmaps_['verticalvelocity']
            norm=colors.SymLogNorm(0.25)
            #contours=np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0]))
            
            plt.contourf(X,Y, Wr, cmap=cmap,norm=norm)
            
            plotting.add_map_locations(extentname, hide_text=True)
            
            plt.colorbar(format=ticker.ScalarFormatter(), pad=0)
            plt.xticks([],[])
            plt.yticks([],[])
            if ii==0:
                plt.title('vertical motion (m/s)')
        
        plt.suptitle("%s weather "%model_version + dtime.strftime("%Y %b %d %H:%M (UTC)"))
        dstamp = dtime.strftime("%Y%m%d%H%M")
        pname="figures/%s/weather_summary/%s/fig_%s%s"%(extentname,model_version,dstamp,ext)
        fio.save_fig(pname,plt)

def read_and_plot_model_run(model_version='waroona_oldold'):
    '''
    '''
    extentname = model_version.split('_')[0]
    extent = plotting._extents_[extentname]
    # font sizes etc
    plotting.init_plots()
    
    if model_version=='waroona_oldold':
        cubes = fio.read_waroona_oldold(extent=extent, add_winds=True)
        u,v,s = cubes.extract(['u','v','s'])
        w, = cubes.extract('upward_air_velocity')
        height = u.coord('Hybrid height').points

        plot_weather_summary(u, v, w, height, 
                             timedim_name='t',
                             model_version=model_version)
    
    elif model_version=='waroona_run1':
        for dtime in fio.model_outputs['waroona_run1']['filedates']:
            cubelists = fio.read_waroona(dtime,extent=extent, add_winds=True)
            u,v,s = cubelists[1].extract(['u','v','s'])
            w, = cubelists[2].extract('upward_air_velocity')
            height = w.coord('level_height').points
            
            plot_weather_summary(u, v, w, height, 
                                 timedim_name='time',
                                 model_version=model_version)
            
#read_and_plot_oldold_run()
#read_and_plot_model_run('waroona_run1')
read_and_plot_model_run('waroona_oldold')
