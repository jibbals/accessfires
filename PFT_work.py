#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 21:04:13 2019

    Try to get fortran PFT calc going
@author: jesse
"""

import numpy as np
from datetime import datetime
import iris
from timeit import default_timer as timer
from metpy.plots import SkewT

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
    print(zsfc.shape, psfc.shape, Tsfc.shape, has_time_dim, len(cubedtimes), latlon)
    print(wwcube)
    print(wwcube.coord('time'))
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

def model_run_PFT_map(mr = 'waroona_run1', hour=datetime(2016,1,5,15)):
    """
    Read model outputs and produce PFT maps for each available time step
    """
    
    plotting.init_plots()
    extentname=mr.split("_")[0]
    extent=plotting._extents_[extentname]
    
    # Read the cubes
    cubes = fio.read_model_run(mr, fdtime=hour,
                               extent=extent,
                               add_z=True, add_RH=True,
                               add_topog=True, add_winds=True,
                               add_theta=True)
    dtimes = utils.dates_from_iris(cubes.extract('u')[0])
    # Read PFT over time,lats,lons (somewhat deresolved..)
    Hskip = 10
    PFT = PFT_from_cubelist(cubes, latskip=Hskip, lonskip=Hskip)
    lats = cubes[0].coord('latitude').points[::Hskip]
    lons = cubes[0].coord('longitude').points[::Hskip]
    terrain = cubes.extract('surface_altitude',strict=True)[::Hskip,::Hskip]
    #print(terrain)
    
    vmin=0
    vmax=10000
    norm = colors.SymLogNorm(200,
                             #vmin=vmin,vmax=vmax
                             ) # linear up to 100, then logarithmic
    for i,dtime in enumerate(dtimes):
        #[print(np.shape(arr)) for arr in [lats,lons,PFT]]
        cs = plt.contourf(lons,lats,PFT[i], 200, 
                          norm=norm, 
                          #vmin=vmin, vmax=vmax, 
                          cmap="jet_r")
        cb = plt.colorbar(cs, 
                          #norm=norm, pad=.015, 
                          format=ticker.ScalarFormatter(),
                          )
        
        ticks = [np.nanmin(PFT[i]), np.nanmax(PFT[i]),0,50,75,100,125,150,175,200,300,400,500,1000,5000,10000]
        ticks.sort()
        cb.set_ticks(ticks)
        #cb.set_clim(0,10000)
        
        ## Add contour lines for terrain height
        cs2 = plt.contour(lons, lats, terrain.data, cmap='terrain', levels=np.arange(0,1001,50))
        plt.clabel(cs2, cs2.levels[0:2], fontsize=13, inline=1, fmt = '%4.0fm') # label lowest two levels
        
        plotting.map_add_locations_extent(extentname,hide_text=True)
        plt.title(dtime.strftime("PFT at %b %d %H:%M (UTC)"))
        fio.save_fig(mr,_sn_,dtime,plt)




if __name__ == '__main__':

    # Run a few tests:
    model_runs=['sirivan_run1','waroona_run1','waroona_old']
    testing = False
    
    for mr in model_runs:
        dtimes = fio.model_outputs[mr]['filedates']
        if testing:
            dtimes = dtimes[0:2]
        for dtime in dtimes:
            model_run_PFT_map(mr=mr,hour=dtime)
    