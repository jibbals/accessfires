#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 21:04:13 2019

    Try to get fortran PFT calc going
@author: jesse
"""

import matplotlib
# don't plot on screen, send straight to file
# this is for lack of NCI display
matplotlib.use('Agg',warn=False)

import numpy as np
from datetime import datetime, timedelta
import iris

import warnings
#import cartopy

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors, ticker

from utilities import fio, utils, plotting


###
## GLOBALS
###
_sn_ = 'PFT_work'


PFT = {'waroona_run2':{'data':None, # manually calculated PFT
                       'units':'Gigawatts',
                       'time':None, # calculated at these times in UTC
                       'latlon':None, # ~ 1km from fire
                       'latlon_stamp':None,
                       'style':'--',
                       'color':'k',
                       },
        'waroona_run1':{'data':np.array([56.5, 42.1, 34.6, 332.6, 53.2]), # manually calculated PFT
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
for mr in fio.model_outputs.keys():
    if mr not in PFT.keys():
        PFT[mr] = PFT['waroona_run2']
        if 'waroona' in mr:
            PFT[mr]['latlon'] = plotting._latlons_['fire_waroona_upwind']
        else:
            PFT[mr]['latlon'] = plotting._latlons_['fire_sirivan_upwind']

## Manual Calculations for PFT:
## YYYYMMDDhhmmrun: 
## 201601060600OLD: q_ML=8.5, th_ML=34, z_fc=680mbar=3.5km, dth=6, U=8m/s :: PFT = .3*3.5**2*6*8 = 176.4 GW
## 201601060800OLD: q_ML=8, th_ML=35, z_fc=685mbar=3.4km, dth=4, U=10.5 :: PFT = .3*3.4**2*4*10.5 = 145.7 GW
## 201601060810RUN1: q_ML=8, th_ML=35.5, z_fc=690mbar=3.33km, dth=2, U=8 :: PFT = .3*3.33**2*2*8 = 53.2 GW
## 201601051500OLD: q_ML=8  , th_ML=35C, star_SP=(10  , 55C), z_fc=690mbar=3.24km, dth=1.5, U=13m/s ::: PFT = .3*3.24**2*1.5*13 = 61.4GW
## 201601060501OLD: q_ML=9  , th_ML=34C, star_SP=(11  , 54C), z_fc=710mbar=3.07km, dth=4.0, U= 6m/s ::: PFT = .3*3.07**2*4*6 = 67.9GW
## 201601051510NEW: q_ML=8  , th_ML=35C, star_SP=(10  , 55C), z_fc=710mbar=3.07km, dth=2.0, U=10m/s ::: PFT = .3*3.07**2*2*10 = 56.5GW

## 201601061210NEW: q_ML=7.2, th_ML=35.5, z_fc = 695mbar=4km, dth=4.5, U=15.4m/s ::: PFt = .3*4**2*4.5*15.4 = 332.6GW

def pft_altitude_vs_pressure(model_run='waroona_run1', latlon=plotting._latlons_['fire_waroona_upwind'],
                             mbar_to_watch=700, datetimes=[datetime(2016,1,5,15)]):
    """
    Retrieve altitudes around 800-600 mbar over a specific lat lon over time.
    This is to calculate z_fc for manual PFT calculation
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

def firepower_comparison(runs=['waroona_run1','waroona_old','waroona_run2','waroona_run3'], localtime=False):
    """
    Plot overlaid time series of two model runs fire power from integral of intensity
    """
    
    lat,lon = PFT[runs[0]]['latlon']
    fpextent = plotting._extents_[runs[0].split('_')[0]]
    
    offset = timedelta(hours=8)
    if lon > 130:
        offset = timedelta(hours=10)
    labels = runs
    #labels = ['new','orig'] # just for poster
    
    nruns=len(runs)
    if nruns>1:
        fig, axes = plt.subplots(len(runs), 1, figsize=[10,8],
                                 sharex=True, sharey=True)
    else:
        fig = plt.figure(figsize=[10,5])
        axes = [plt.gca(),]
    
    for i,run in enumerate(runs):
        ax=axes[i]
        plt.sca(ax)
        # read all the fire data
        sh, = fio.read_fire(model_run=run, dtimes=None,
                            extent=fpextent, firefront=False, sensibleheat=True)
        assert sh is not None, "Missing sensible heat file"
        
        ftimes = utils.dates_from_iris(sh)
        if localtime:
            ftimes = np.array([ft + offset for ft in ftimes ])
        firepower = np.sum(utils.firepower_from_cube(sh), axis=(1,2)) # over time
        
        # remove zeros:
        prefire = np.isclose(np.cumsum(firepower), 0)
        firepower[prefire] = np.NaN
        
        ## Plot firepower
        plt.plot_date(ftimes, firepower, '-', label='firepower')
        
        ## Read PFT
        # calculated using kevin's code
        pft, ptimes, _, _ = fio.read_pft(run,lats=lat,lons=lon)
        if localtime:
            ptimes = np.array([pt + offset for pt in ptimes ])
        
        ## plot PFT
        plt.plot_date(ptimes, pft, '--', label='PFT')
        
        ## Fix xticks and stuff for multiple plots
        if i < nruns-1:
            # remove x ticks and labels,
            plt.xticks([],[])
        if i == 0:
            plt.title("Firepower and PFT")
            plt.legend(loc='best')
        plt.title(labels[i],y=0.9)
        
    # fix up labels, axes
    if nruns>1:
        fig.subplots_adjust(hspace=0.0)
    
    # format x-ticks date times
    ax = plt.gca()
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
    ax.xaxis.set_minor_locator(mdates.HourLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    fio.save_fig(run, _sn_, 'firepower.png',subdir='comparison',plt=plt)

def PFT_map(PFT,plats,plons, vmax=250):
    """
    Plot top-down map of PFT with contour and optional colour bar
    vmax - 0 is linear reds
    vmax+ is shown as light blue
    """
    
    # remove negatives before plotting log scale (or else warnings apear)
    PFT_pos = np.copy(PFT)
    #PFT_pos[PFT_pos > vmax] = vmax
    PFT_pos[PFT_pos<=0] = np.NaN
    
    cs = plt.contourf(plons, plats, PFT_pos,
                      levels = np.arange(0,vmax,30),
                      cmap="Reds_r",
                      extend="max",
                      vmax=vmax, vmin=0)
    cs.cmap.set_over("blue")
    cs.changed()
    
    tickvals = np.append(np.arange(0,vmax,100),vmax)
    #print("DEBUG:",tickvals)
    cb = plt.colorbar(cs,pad=0.01,ticks=tickvals)
    cb.set_label('PFT [Gigawatts]')
    
    cb.set_ticks(list(tickvals))
    cb.ax.set_yticklabels(list(tickvals))
    
    
    return cs, cb
    
def model_run_PFT_summary(model_run='waroona_run1', hour=datetime(2016,1,5,15)):
    '''
    Show PFT map, underneath topography overlaid with curly surface wind map
    '''
    extentname=model_run.split('_')[0]
    extent = plotting._extents_[extentname]
    
    cubes = fio.read_model_run(model_run,[hour],extent=extent, 
                               add_topog=True,add_winds=True)
    u0,v0 = cubes.extract(['u','v'], strict=True)
    dtimes = utils.dates_from_iris(u0,remove_seconds=True)
    lats = u0.coord('latitude').points
    lons = u0.coord('longitude').points
    heights = utils.height_from_iris(u0)
    surface = heights < 500
    u=np.mean(u0[:,surface,:,:].data, axis=1)
    v=np.mean(v0[:,surface,:,:].data, axis=1)
    
    terrain0 = cubes.extract('surface_altitude',strict=True)
    terrain = terrain0.data
    
    ## Read PFT
    pft, ptimes, plats, plons = fio.read_pft(model_run, dtimes, lats, lons)
    
    ## Read fire front
    ff, = fio.read_fire(model_run, dtimes, extent=extent)
    
    
    for i,dtime in enumerate(dtimes):
        ## First subplot will be topography and locations
        plt.subplot(2,1,1)
        plotting.map_topography(extent, terrain, lats, lons, title='')
        plotting.map_add_locations_extent(extentname,hide_text=False)
        # add fire front
        if ff is not None:
            plotting.map_fire(ff[i].data,lats,lons)
        
        ## Second plot will be PFT map with 'surface' winds overlaid
        plt.subplot(2,1,2)
        PFT_map(pft[i],plats,plons)
        
        # overlay winds
        plt.streamplot(lons,lats,u[i],v[i], color='k', 
                       density=(.8, .5))
        
        plotting.map_add_locations_extent(extentname,hide_text=True)
        ## Turn off the tick values
        plt.xticks([]); plt.yticks([])
        ## Title, save, close
        plt.suptitle(dtime.strftime("PFT and mean winds below 500m at %b %d %H:%M (UTC)"))
        fio.save_fig(model_run,_sn_,dtime,plt)


if __name__ == '__main__':
    
    ## Compare firepower/PFT for some runs
    #firepower_comparison(runs=['waroona_old','waroona_run1'])
    #firepower_comparison()

    ## Summary figure for PFT at a site for one output hour
    if False:
        for mr in ['sirivan_run1','sirivan_run1_hr']:
            #['waroona_run3', 'waroona_run1']:
            dtimes = fio.model_outputs[mr]['filedates']
            for hour in np.arange(5,20):
                model_run_PFT_summary(model_run=mr, hour=dtimes[hour])
    
    
    if True:
        firepower_comparison(['sirivan_run1',])
    
