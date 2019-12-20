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
from datetime import datetime
import iris
from timeit import default_timer as timer
from metpy.plots import SkewT
import warnings
import cartopy

import matplotlib.pyplot as plt
from matplotlib import colors, ticker

from utilities import fio, utils, plotting

## PFT is my port of Kevin's proprietry fortran code
## if this is not available, we are out of luck
import utilities.fortran.PFT as HFj    


###
## GLOBALS
###
_sn_ = 'PFT_work'



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

def PFT_map(PFT,plats,plons,colorbar=True, lines=[100]):
    """
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
    
def model_run_PFT_summary(mr='waroona_run1', hour=datetime(2016,1,5,15)):
    '''
    Show PFT map, underneath topography overlaid with curly surface wind map
    '''
    extentname=mr.split('_')[0]
    extent = plotting._extents_[extentname]
    
    cubes = fio.read_model_run(mr,[hour],extent=extent, 
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
    pft, ptimes, plats, plons = fio.read_pft(mr, dtimes, lats, lons)
    
    ## Read fire front
    hasfires= mr in ['waroona_run1','waroona_old','waroona_run2','sirivan_run1']
    if hasfires:
        ff, = fio.read_fire(mr, dtimes, extent=extent)
    
    
    for i,dtime in enumerate(dtimes):
        ## First subplot will be topography and locations
        plt.subplot(2,1,1)
        plotting.map_topography(extent, terrain, lats, lons, title='')
        plotting.map_add_locations_extent(extentname,hide_text=False)
        # add fire front
        if hasfires:
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
        fio.save_fig(mr,_sn_,dtime,plt)


if __name__ == '__main__':
    
    ## test method
    #model_run_PFT_summary()
    
    run_everything=True
    if run_everything:
        for mr in ['waroona_run2','waroona_run2uc']:#['sirivan_run1','waroona_run1','waroona_old']:
            dtimes = fio.model_outputs[mr]['filedates']
            for dtime in dtimes:
                model_run_PFT_summary(mr=mr, hour=dtime)
    
