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
from timeit import default_timer as timer
import warnings
#import cartopy

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import colors, ticker

from utilities import fio, utils, plotting

## PFT is my port of Kevin's proprietry fortran code
## if this is not available, we are out of luck
import utilities.fortran.PFT as HFj    


###
## GLOBALS
###
_sn_ = 'PFT_work'


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

## Manual Calculations for PFT:
## YYYYMMDDhhmmrun: 
## 201601060600OLD: q_ML=8.5, th_ML=34, z_fc=680mbar=3.5km, dth=6, U=8m/s :: PFT = .3*3.5**2*6*8 = 176.4 GW
## 201601060800OLD: q_ML=8, th_ML=35, z_fc=685mbar=3.4km, dth=4, U=10.5 :: PFT = .3*3.4**2*4*10.5 = 145.7 GW
## 201601060810RUN1: q_ML=8, th_ML=35.5, z_fc=690mbar=3.33km, dth=2, U=8 :: PFT = .3*3.33**2*2*8 = 53.2 GW
## 201601051500OLD: q_ML=8  , th_ML=35C, star_SP=(10  , 55C), z_fc=690mbar=3.24km, dth=1.5, U=13m/s ::: PFT = .3*3.24**2*1.5*13 = 61.4GW
## 201601060501OLD: q_ML=9  , th_ML=34C, star_SP=(11  , 54C), z_fc=710mbar=3.07km, dth=4.0, U= 6m/s ::: PFT = .3*3.07**2*4*6 = 67.9GW
## 201601051510NEW: q_ML=8  , th_ML=35C, star_SP=(10  , 55C), z_fc=710mbar=3.07km, dth=2.0, U=10m/s ::: PFT = .3*3.07**2*2*10 = 56.5GW

## 201601061210NEW: q_ML=7.2, th_ML=35.5, z_fc = 695mbar=4km, dth=4.5, U=15.4m/s ::: PFt = .3*4**2*4.5*15.4 = 332.6GW

def PFT_from_cubelist(cubes0, latlon=None, tskip=None, latskip=None, lonskip=None):
    """
    Wrapper to call the PFT function using cubelist of data
    
    Loop over required dimensions, 200 PFT calculations take ~ 1 minutes.
    Subsetting or striding dimensions is recommended for larger data sets.
    cubelist requires these named cubes:
        'air_temperature', 'specific_humidity', 'potential_temperature',
        'air_pressure','u','v','upward_air_velocity',
        'surface_altitude', 'surface_air_pressure', 'surface_temperature'
    if latlon is specified, just one point is extracted using linear interpolation
    
    """
    # local cubelist copy, so that original cubes aren't modified
    cubes = iris.cube.CubeList()
    # if there's no surface data, use first level as surf
    if len(cubes0.extract('surface_temperature')) == 0:
        st = cubes0.extract('air_temperature')[0].copy()
        st.rename("surface_temperature")
        st = st[:,0,:,:]
        cubes0.append(st)
    if len(cubes0.extract('surface_air_pressure')) == 0:
        sap = cubes0.extract('air_pressure')[0].copy()
        sap.rename("surface_air_pressure")
        sap = sap[:,0,:,:]
        cubes0.append(sap)
        
    for cube in cubes0.extract(['air_temperature', 'specific_humidity', 
                    'potential_temperature', 'air_pressure', 'u', 'v',
                    'upward_air_velocity', 'surface_altitude', 
                    'surface_air_pressure', 'surface_temperature'], strict=True):
        cubes.append(cube.copy())
    
    
    ## first interpolate everything to latlon
    if latlon is not None:
        lat,lon = latlon
        for i in range(len(cubes)):
            cubes[i] = cubes[i].interpolate([('longitude',lon), ('latitude',lat)],
                                            iris.analysis.Linear())
    
    has_time_dim = len(cubes.extract('u')[0].coords('time')) == 1
    
    # Now subset everything based on skips
    if ((latskip is not None) or (lonskip is not None)) and (latlon is None):
        latslice = slice(None, None, latskip)
        lonslice = slice(None, None, lonskip)
        
        for i in range(len(cubes)):
            # cube needs to have spacial grid
            if len(cubes[i].coords('longitude')) == 1:
                # cube may have vertical grid and or temporal grid
                cshape=cubes[i].shape
                if len(cshape) == 2:
                    cubes[i] = cubes[i][latslice,lonslice]
                elif len(cshape) == 3:
                    cubes[i] = cubes[i][:,latslice,lonslice]
                elif len(cshape) == 4:
                    cubes[i] = cubes[i][:,:,latslice,lonslice]
    if has_time_dim and (tskip is not None):
        tslice = slice(None, None, tskip)
        # for each cube
        for i in range(len(cubes)):
            # if the cube has a time dim, of length greater than 1, slice it
            if len(cubes[i].coords('time')) == 1:
                if len(cubes[i].coord('time').points) > 1:
                    cubes[i] = cubes[i][tslice]
    
    if has_time_dim:
        cubedtimes = utils.dates_from_iris(cubes.extract('u')[0])
    
    # now for easy reading pull out cubes
    TTcube, qqcube, thcube = cubes.extract(['air_temperature', 
                                            'specific_humidity',
                                            'potential_temperature'],strict=True)
    prcube, uucube, vvcube = cubes.extract(['air_pressure', 'u','v'],strict=True)
    wwcube = cubes.extract('upward_air_velocity',strict=True)
    
    # surface metrics
    # surface values in old run are on different time dimension...!?!
    zsfc, psfc, Tsfc = cubes.extract(['surface_altitude', 
                                      'surface_air_pressure', 
                                      'surface_temperature'])
    #print(zsfc.shape, psfc.shape, Tsfc.shape, has_time_dim, len(cubedtimes), latlon)
    #print(wwcube)
    #print(wwcube.coord('time'))
    zsfc = zsfc.data # m
    if len(zsfc.shape) == 0:
        zsfc = float(zsfc)
    
    # psfc and Tsfc may not have time dim.
    # if they do, or if nothing has time dims, just treat normally
    #if (len(wwcube.coord('time').points) == 1) or (not has_time_dim):
    #    psfc = psfc.data # Pa
    #    Tsfc = Tsfc.data # K
    #    if len(psfc.shape) == 0:
    #        psfc = float(psfc)
    #        Tsfc = float(Tsfc)
    # if they don't, and our other data has a time dim, repeat these over the time dim
    if psfc.shape[0] != wwcube.shape[0]:
        # repeat along time dim
        if latlon is None:
            psfc = np.repeat(psfc.data[np.newaxis,:,:], len(cubedtimes),axis=0)
            Tsfc = np.repeat(Tsfc.data[np.newaxis,:,:], len(cubedtimes),axis=0)
        else:
            psfc = np.repeat(float(psfc.data), len(cubedtimes),axis=0)
            Tsfc = np.repeat(float(Tsfc.data), len(cubedtimes),axis=0)
    else:
        psfc = psfc.data # Pa
        Tsfc = Tsfc.data # K
        if len(psfc.shape) == 0:
            psfc = float(psfc)
            Tsfc = float(Tsfc)
            
    # Return array
    PFT = np.zeros(Tsfc.shape)+np.NaN 
    
    # let's time how long it takes
    start = timer()
    
    ## Loop over time, lat, lon dimensions
    if latlon is None:
        lats,lons = TTcube.coord('latitude').points,TTcube.coord('longitude').points
        for yi in range(len(lats)):
            for xi in range(len(lons)):
                zsfc0 = zsfc[yi,xi] # zsfc never has time dim
                if has_time_dim:
                    for ti in range (len(cubedtimes)):
                        TT, qq = TTcube[ti,:,yi,xi].data.data, qqcube[ti,:,yi,xi].data.data
                        uu,vv,ww = uucube[ti,:,yi,xi].data.data, vvcube[ti,:,yi,xi].data.data, wwcube[ti,:,yi,xi].data.data
                        th,pr = thcube[ti,:,yi,xi].data.data, prcube[ti,:,yi,xi].data.data
                        psfc0, Tsfc0 = psfc[ti,yi,xi], Tsfc[ti,yi,xi]

                        # Find pressure level at which T = Tmin (-20 degr C)
                        # get first instance where TT is less than Tmin
                        Tmin_indices = TT < 253.15
                        pmin = pr[Tmin_indices][0]
                        
                        frets = HFj.PFT(TT,qq,uu,vv,ww,th,pr, 
                                        zsfc0,psfc0,Tsfc0,
                                        Pmin=pmin)
                        PFT[ti,yi,xi] = frets[8]/1e9 # G Watts
                # if no time dimension
                else:
                    TT, qq = TTcube[:,yi,xi].data.data, qqcube[:,yi,xi].data.data
                    uu,vv,ww = uucube[:,yi,xi].data.data, vvcube[:,yi,xi].data.data, wwcube[:,yi,xi].data.data
                    th,pr = thcube[:,yi,xi].data.data, prcube[:,yi,xi].data.data
                    psfc0, Tsfc0 = psfc[yi,xi], Tsfc[yi,xi]

                    # Find pressure level at which T = Tmin (-20 degr C)
                    # get first instance where TT is less than Tmin
                    Tmin_indices = TT < 253.15
                    pmin = pr[Tmin_indices][0]
                    
                    frets = HFj.PFT(TT,qq,uu,vv,ww,th,pr, 
                                    zsfc0,psfc0,Tsfc0,
                                    Pmin=pmin)
                    PFT[yi,xi] = frets[8]/1e9 # G Watts
    # If we have no lat, lon dim
    else:
        zsfc0 = zsfc
        if has_time_dim:
            for ti in range (len(cubedtimes)):
                TT, qq = TTcube[ti,:].data.data, qqcube[ti,:].data.data
                uu,vv,ww = uucube[ti,:].data.data, vvcube[ti,:].data.data, wwcube[ti,:].data.data
                th,pr = thcube[ti,:].data.data, prcube[ti,:].data.data
                psfc0, Tsfc0 = psfc[ti], Tsfc[ti]

                # Find pressure level at which T = Tmin (-20 degr C)
                # get first instance where TT is less than Tmin
                Tmin_indices = TT < 253.15
                pmin = pr[Tmin_indices][0]
                
                frets = HFj.PFT(TT,qq,uu,vv,ww,th,pr, 
                                zsfc0,psfc0,Tsfc0,
                                Pmin=pmin)
                PFT[ti] = frets[8]/1e9 # G Watts
        # if no time dimension and no spatial dim
        else:
            TT, qq = TTcube.data.data, qqcube.data.data
            uu,vv,ww = uucube.data.data, vvcube.data.data, wwcube.data.data
            th,pr = thcube.data.data, prcube.data.data
            psfc0, Tsfc0 = psfc, Tsfc

            # Find pressure level at which T = Tmin (-20 degr C)
            # get first instance where TT is less than Tmin
            Tmin_indices = TT < 253.15
            pmin = pr[Tmin_indices][0]
            
            frets = HFj.PFT(TT,qq,uu,vv,ww,th,pr, 
                            zsfc0,psfc0,Tsfc0,
                            Pmin=pmin)
            PFT = frets[8]/1e9 # G Watts
    end = timer()
    print("Info: time to produce PFT(%s): %.2f minutes"%(str(PFT.shape), (end-start)/60.0))
    return PFT

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

def firepower_comparison(runs=['waroona_run1','waroona_old','waroona_run2'], localtime=False):
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

def PFT_map(PFT,plats,plons,colorbar=True, lines=[100]):
    """
    Plot top-down map of PFT with contour and optional colour bar
    """
    cnorm = colors.SymLogNorm(1,vmin=0,vmax=1000)
    levs = np.union1d([0],np.logspace(0,3,20))
    
    cs = plt.contourf(plons,plats,PFT, levs, 
                      norm=cnorm,
                      locator=ticker.LogLocator(),
                      cmap="YlOrRd_r")
    cb = None
    if colorbar:
        cb = plt.colorbar(pad=0.01)
        cb.set_label('PFT [Gigawatts]')
        cb.set_ticks([1e1, 1e2, 1e3])
        cb.ax.set_yticklabels(['10$^1$','10$^2$','10$^3$'])
    
    if lines is not None:
        with warnings.catch_warnings():
            # ignore warning when there are no fires:
            warnings.simplefilter('ignore')
            plt.contour(plons,plats,PFT,np.array(lines))
        
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
    heights = u0.coord('level_height').points
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
            with warnings.catch_warnings():
                # ignore warning when there are no fires:
                warnings.simplefilter('ignore')
                plt.contour(lons,lats,np.transpose(ff[i].data),np.array([0]), colors='red')
        
        # add scale
        #plotting.scale_bar(plt.gca(), cartopy.crs.PlateCarree(),10)
        
        ## Second plot will be PFT map with 'surface' winds overlaid
        plt.subplot(2,1,2)
        PFT_map(pft[i],plats,plons)
        # overlay winds
        plotting.map_quiver(u[i].data,v[i].data,lats,lons, scale=140)
        plotting.map_add_locations_extent(extentname,hide_text=True)
        ## Turn off the tick values
        plt.xticks([]); plt.yticks([])
        ## Title, save, close
        plt.suptitle(dtime.strftime("PFT and mean winds below 500m at %b %d %H:%M (UTC)"))
        fio.save_fig(model_run,_sn_,dtime,plt)


if __name__ == '__main__':
    
    ## test method
    #model_run_PFT_summary()
    
    firepower_comparison(runs=['waroona_old','waroona_run1'])
    
    
    #run_everything=True
    #if run_everything:
    #    for mr in ['waroona_run2','waroona_run2uc']:#['sirivan_run1','waroona_run1','waroona_old']:
    #        dtimes = fio.model_outputs[mr]['filedates']
    #        for dtime in dtimes:
    #            model_run_PFT_summary(model_run=mr, hour=dtime)
    
