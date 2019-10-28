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

from utilities import fio, utils, plotting


## PFT is my port of Kevin's proprietry fortran code
## if this is not available, we are out of luck
import utilities.fortran.PFT as HFj    

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
    for cube in cubes0.extract(['air_temperature', 'specific_humidity', 
                    'potential_temperature', 'air_pressure', 'u', 'v',
                    'upward_air_velocity', 'surface_altitude', 
                    'surface_air_pressure', 'surface_temperature'], strict=True):
        cubes.append(cube.copy())
    
    ## first interpolate everything to latlon
    if latlon is not None:
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
    zsfc = zsfc.data # m
    if len(zsfc.shape) == 0:
        zsfc = float(zsfc)
    
    # psfc and Tsfc may not have time dim.
    # if they do, or if nothing has time dims, just treat normally
    if (len(psfc.coords('time')) == 1) or (not has_time_dim):
        psfc = psfc.data # Pa
        Tsfc = Tsfc.data # K
        if len(psfc.shape) == 0:
            psfc = float(psfc)
            Tsfc = float(Tsfc)
    # if they don't, and our other data has a time dim, repeat these over the time dim
    else:
        # repeat along time dim
        if latlon is None:
            psfc = np.repeat(psfc.data[np.newaxis,:,:], len(cubedtimes),axis=0)
            Tsfc = np.repeat(Tsfc.data[np.newaxis,:,:], len(cubedtimes),axis=0)
        else:
            psfc = np.repeat(float(psfc.data), len(cubedtimes),axis=0)
            Tsfc = np.repeat(float(Tsfc.data), len(cubedtimes),axis=0)
    
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

if __name__ == '__main__':

    # Run a few tests:
    
    ## First read a profile, and surface information required by fortran subroutine
    latlon = plotting._latlons_['fire_waroona_upwind']
    lat,lon = latlon
    extent = [lon-.01, lon+.01, lat-.01, lat+.01] # just grab real close to latlon
    cubes=fio.read_model_run('waroona_run1',fdtime=[datetime(2016,1,6,5),datetime(2016,1,6,6)],
                             extent=extent, add_dewpoint=True, 
                             add_winds=True, add_RH=False,
                             add_topog=True, add_z=True, add_theta=True)
    
    PFT0 = PFT_from_cubelist(cubes) # whole map with time steps etc
    
    PFT1 = PFT_from_cubelist(cubes,latlon=latlon) # just a time series
    
    PFT2 = PFT_from_cubelist(cubes,latskip=2,tskip=2) # subset
    
    print("OUTPUT:           Min,    mean,    max for PFT calculations")
    print("OUTPUT: PFT0:  %.2f, %.2f, %.2f"%(np.min(PFT0), np.mean(PFT0), np.max(PFT0)))
    print("OUTPUT: PFT1:  %.2f, %.2f, %.2f"%(np.min(PFT1), np.mean(PFT1), np.max(PFT1)))
    print("OUTPUT: PFT2:  %.2f, %.2f, %.2f"%(np.min(PFT2), np.mean(PFT2), np.max(PFT2)))
    plt.subplot(211)
    for yi in range(PFT0.shape[1]):
        for xi in range(PFT0.shape[2]):
            plt.plot(PFT0[:,yi,xi], alpha=0.4, linestyle=':')
    for yi in range(PFT2.shape[1]):
        for xi in range(PFT2.shape[2]):
            plt.plot(PFT2[:,yi,xi], alpha=0.5, linestyle='--')

    plt.plot(PFT1, linewidth=2, color='k') # should be a series
    plt.title("PFT series at, or near, %.3f, %.3f"%(lat,lon))
    
    plt.subplot(212)
    plt.contourf(PFT0[0])
    plt.colorbar()