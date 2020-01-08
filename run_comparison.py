#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 15:45:40 2019
    
    Compare two model runs together...
    
@author: jesse
"""

import matplotlib
#matplotlib.use("Agg",warn=False)

from matplotlib import colors, ticker

import numpy as np
import matplotlib.pyplot as plt
import warnings
from datetime import datetime
from scipy.stats import gaussian_kde, cumfreq

from utilities import utils, plotting, fio, constants

## GLOBAL
#Script name
_sn_ = 'run_comparison'

def compare_surface(mr1='waroona_run2', mr2='waroona_run2uc', hour=datetime(2016,1,5,15)):
    """
        3 by 3 plot, ROWS: H winds, V winds, Theta, COLS: Run2, Run2UC, DIFF
    """
    fig, axes = plt.subplots(3,3)
    

def compare_winds(mr1='waroona_run2', mr2='waroona_run2uc', hour=datetime(2016,1,5,15)):
    """
        Compare winds between two different runs, looking at vertically binned averages
        first row shows model run 1 (and figures are saved under this run)
        second row shows model run 2
        final row(s) show density or cumulative density from averaged metric
        
        two separate folders:
            horizontal/
                4 averaged binds (altitudes 0-500, 500-2000, 2000-5000, 5000-9000)
                mr1
                mr2
                distribution histograms for wind speed and direction
            vertical
                same as horizontal, with cdf for wind speed
    """
    extentname = mr1.split('_')[0]
    extent = plotting._extents_[extentname]
    
    ## Read a model run
    cubes1 = fio.read_model_run(mr1, hour, extent=extent, add_winds=True)
    cubes2 = fio.read_model_run(mr2, hour, extent=extent, add_winds=True)
    
    # pull out horizontal winds
    u,v,s,wd = cubes1.extract(['u','v','s','wind_direction'])
    cu,cv,cs,cwd = cubes2.extract(['u','v','s','wind_direction'])
    dates = utils.dates_from_iris(u)
    height = s.coord('level_height').points
    lats = s.coord('latitude').points
    lons = s.coord('longitude').points
    ff1,ff2=None,None
    if fio.model_outputs[mr1]['hasfire']:
        ff1, = fio.read_fire(mr1,dtimes=dates,extent=extent,firefront=True)
    if fio.model_outputs[mr2]['hasfire']:
        ff2, = fio.read_fire(mr2,dtimes=dates,extent=extent,firefront=True)
    
    # density plot arguments
    bandwidth=1
    # x axis for wind speed
    xs = np.linspace(0,20,100)
    # x axis for wind dir
    xwd = np.linspace(0,360,100)
    # x axis for vert wind speed
    xw = np.linspace(-32,32,100)
    
    ## Colourmap setup for contourf plots
    # Horizontal wind colorbar should be the same for both runs
    hmaxthresh=[5,10,15,20,25,30,35]
    hcmap=plotting._cmaps_['windspeed']
    # vertical wind colourbar is constant
    wcmap=plotting._cmaps_['verticalvelocity']
    wnorm=colors.SymLogNorm(0.25) # linear to +- 0.25, then log scale
    wcontours=np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0]))
    
    # make 4 vertical bins
    row1 = (0<=height) * (height<500)
    row2 = (500<=height) * (height<2000)
    row3 = (2000<=height) * (height<5000)
    row4 = (5000<=height) * (height<9000)
    
    # loop over datetimes:
    for di, date in enumerate(dates):
        
        s1, s2, s3, s4 = [np.mean(s[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        u1, u2, u3, u4 = [np.mean(u[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        v1, v2, v3, v4 = [np.mean(v[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        wd1, wd2, wd3, wd4 = [np.mean(wd[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        
        cs1, cs2, cs3, cs4 = [np.mean(cs[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        cu1, cu2, cu3, cu4 = [np.mean(cu[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        cv1, cv2, cv3, cv4 = [np.mean(cv[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        cwd1, cwd2, cwd3, cwd4 = [np.mean(cwd[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        
        
        # Plotting
        plt.close()
        fig, axes = plt.subplots(4,4,figsize=[12,12])
        for i, (si, ui, vi, wdi, csi, cui, cvi, cwdi) in \
            enumerate(zip([s1,s2,s3,s4],[u1,u2,u3,u4],[v1,v2,v3,v4],[wd1,wd2,wd3,wd4],
                  [cs1,cs2,cs3,cs4],[cu1,cu2,cu3,cu4],[cv1,cv2,cv3,cv4],[cwd1,cwd2,cwd3,cwd4])):
            ## Show contourf of wind speed
            plt.sca(axes[0,i])
            
            # first determine colourmap contours so that they are the same between plots (and useful)
            hmax_index=np.sum(np.max([si,csi])>hmaxthresh)
            hcontours=np.linspace(0,hmaxthresh[hmax_index],20)
            
            # plot the filled contour for h-wind speeds
            plotting.map_contourf(extent, si,lats,lons,cmap=hcmap,clabel="",clevs=hcontours,cbar=False,cbarform=None)
            # overlaid with quiver of wind dir
            plotting.map_quiver(ui,vi,lats,lons,nquivers=7)
            # Add locations and fire
            plotting.map_add_locations_extent(extentname, hide_text=True)
            if ff1 is not None:
                plotting.map_fire(ff1[di].data,lats,lons)
            # Add ylabel on left most plot
            if i==0: plt.ylabel(mr1)
            # add title along top row
            plt.title(['<500m','500m-2000m','2km-5km','5km-9km'][i])
            
            ## for comparison model also
            plt.sca(axes[1,i])
            img,_ = plotting.map_contourf(extent, csi, lats, lons,cmap=hcmap,clabel="",clevs=hcontours,cbar=False,cbarform=None)
            
            # overlaid with quiver of wind dir
            plotting.map_quiver(cui,cvi,lats,lons,nquivers=7)
            # add fire
            if ff2 is not None:
                plotting.map_fire(ff2[di].data,lats,lons)
            # add label on leftmost
            if i==0: plt.ylabel(mr2)
            
            ## Add colourbar for column
            cbar_ax = fig.add_axes([0.14+0.2025*i, 0.508, 0.14, 0.01])# X Y Width Height
            cbar=fig.colorbar(img, cax=cbar_ax, format=ticker.ScalarFormatter(), pad=0,
                              orientation='horizontal')
            cbar.set_ticks(np.arange(0,hmaxthresh[hmax_index]+1,5))
            cbar.set_ticklabels(np.arange(0,hmaxthresh[hmax_index]+1,5))
            
            ## add density plot for wind speed
            plt.sca(axes[2,i])
            # density function:
            sdens = gaussian_kde(si.flatten(),bw_method=bandwidth)
            plt.plot(xs,sdens(xs), label=mr1,linewidth=2)
            csdens = gaussian_kde(csi.flatten(),bw_method=bandwidth)
            plt.plot(xs,csdens(xs), label=mr2)
            plt.yticks([],[])
            if i==0: plt.ylabel('wind speed density')
            if i==3: plt.legend()
            
            
            ## add density plot for wind dir
            plt.sca(axes[3,i])
            # density function:
            wddens = gaussian_kde(wdi.flatten(),bw_method=bandwidth)
            plt.plot(xwd,wddens(xwd), label=mr1, linewidth=2)
            cwddens = gaussian_kde(cwdi.flatten(),bw_method=bandwidth)
            plt.plot(xwd,cwddens(xwd), label=mr2, linewidth=2, linestyle='--')
            plt.yticks([],[])
            if i==0: plt.ylabel('wind dir density')
            if i==3: plt.legend()
            
        plt.suptitle(date.strftime("%Y%m%d %H:%M(UTC)"))
        fio.save_fig(mr1, _sn_, date, plt, subdir='horizontal')
    
        ## Also want to look at vertical winds
        w, = cubes1.extract(['upward_air_velocity'])
        cw, = cubes2.extract(['upward_air_velocity'])
        
        w1, w2, w3, w4 = [np.mean(w[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        cw1, cw2, cw3, cw4 = [np.mean(cw[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        
        # Plotting
        plt.close()
        fig, axes = plt.subplots(3,4,figsize=[12,12])
        for i, (wi, cwi) in enumerate(zip([w1,w2,w3,w4],[cw1,cw2,cw3,cw4])):
            ## Show contourf of wind speed
            plt.sca(axes[0,i])
            plotting.map_contourf(extent, wi, lats,lons,cmap=wcmap,clabel="",clevs=wcontours,norm=wnorm,cbar=False,cbarform=None)
            plotting.map_add_locations_extent(extentname, hide_text=True)
            if ff1 is not None:
                plotting.map_fire(ff1[di].data,lats,lons)
            if i==0: plt.ylabel(mr1)
            plt.title(['<500m','500m-2000m','2km-5km','5km-9km'][i])
            
            ## for comparison model also
            plt.sca(axes[1,i])
            img,_=plotting.map_contourf(extent, cwi, lats,lons,cmap=wcmap,clabel="",clevs=wcontours,norm=wnorm,cbar=False,cbarform=None)
            plotting.map_add_locations_extent(extentname, hide_text=True)
            if ff2 is not None:
                plotting.map_fire(ff2[di].data,lats,lons)
            if i==0: plt.ylabel(mr2)
            
            ## add density plot for wind speed
            plt.sca(axes[2,i])
            ## density function:
            #wdens = gaussian_kde(wi.flatten(),bw_method=bandwidth)
            #plt.plot(xw,wdens(xw), label=mr1, linewidth=2)
            #cwdens = gaussian_kde(cwi.flatten(),bw_method=bandwidth)
            #plt.plot(xw,cwdens(xw), label=mr2, linewidth=2, linestyle='--')
            ## cumulative frequency
            wcdf = cumfreq(wi.flatten(), numbins=len(xw), defaultreallimits=(min(xw),max(xw)))
            plt.plot(xw,wcdf.cumcount,label=mr1, linewidth=2)
            cwcdf = cumfreq(cwi.flatten(), numbins=len(xw), defaultreallimits=(min(xw),max(xw)))
            plt.plot(xw,cwcdf.cumcount,label=mr2, linewidth=2, linestyle='--')
            plt.yticks([],[])
            if i==0: plt.ylabel('wind speed cdf')
            if i==3: plt.legend()
        ## Add colourbar
        cbar_ax = fig.add_axes([.35, 0.37, 0.3, 0.015])# X Y Width Height
        cbar=fig.colorbar(img, cax=cbar_ax, format=ticker.ScalarFormatter(), pad=0,
                          orientation='horizontal')
        cbar.set_ticks([-32,-8,-2,0,2,8,32])
        cbar.set_ticklabels([-32,-8,-2,0,2,8,32])
        # title and save
        plt.suptitle(date.strftime("%Y%m%d %H:%M(UTC)"))
        fio.save_fig(mr1, _sn_, date, plt, subdir='vertical')
    

def compare_at_site(mr1='waroona_run2', mr2='waroona_run2uc', latlon = plotting._latlons_['AWS_wagerup']):
    """
    look at time series of metrics from two runs, compared to AWS if possible
    """
    print("TBD")
        



if __name__=='__main__':
    
    #compare_winds()
    for fdate in fio.model_outputs['waroona_run2']['filedates']:
        compare_winds(hour=fdate)
    print("run_comparison.py done")
