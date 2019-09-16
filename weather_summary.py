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
                         ext='.png', timedim_name='time',
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
    
    for ti, dtime in enumerate(dtimes):

        row1 = (100<=height) * (height<500)
        row2 = (500<=height) * (height<1500)
        row3 = (1500<=height) * (height<3000)
        row4 = (3000<=height) * (height<5000)
        
        plt.figure(figsize=[10,10])
        for ii,row in enumerate([row1,row2,row3,row4]):
            plt.subplot(4,2,ii*2+1)
            
            Ui = U[ti,row]
            # only use interpolated data: this is just for oldold run on local machine...
            Vi = V[ti,row].interpolate([('latitude',V.coord('latitude').points)],
                                   iris.analysis.Linear())
            
            # mean of levels
            Ur = np.mean(Ui.data,axis=0)
            Vr = np.mean(Vi.data,axis=0)
            windspeed = np.hypot(Ur,Vr)
            
            # value skip to vectors aren't so dense
            vsu, vsv=30, 50 # lonskip, latskip
            skip = (slice(None,None,vsv),slice(None,None,vsu))
            
            Y=Ui.coord('latitude').points
            X=Vi.coord('longitude').points
            Uvs,Vvs = Ur[skip].data, Vr[skip].data
            
            # Normalize the arrows:
            Uvs = Uvs / np.sqrt(Uvs**2 + Vvs**2);
            Vvs = Vvs / np.sqrt(Uvs**2 + Vvs**2);
            
            plt.contourf(X, Y, windspeed.data)
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
            
            plt.contourf(X,Y,Wr, cmap=cmap,norm=norm)
            plt.colorbar(format=ticker.ScalarFormatter(), pad=0)
            plt.xticks([],[])
            plt.yticks([],[])
            if ii==0:
                plt.title('vertical motion (m/s)')
        
        plt.suptitle("(oldold)  "+dtime.strftime("%Y %b %d %H:%M (UTC)"))
        dstamp = dtime.strftime("%Y%m%d%H%M")
        extentname=['sirivan','waroona'][dtime.year==2016]
        pname="figures/%s/weather_summary/oldold/fig_%s%s"%(extentname,dstamp,ext)
        fio.save_fig(pname,plt)



def read_and_plot_oldold_run(xwind_path = 'data/waroona_oldold/oldold_xwind_s5_subset.nc',
                             ywind_path = 'data/waroona_oldold/oldold_ywind_s5_subset.nc',
                             zwind_path = 'data/waroona_oldold/oldold_zwind_s5_subset.nc'):
    '''
    '''
    extentname = 'waroona'
    extent=plotting._extents_[extentname]
    West,East,South,North = extent
    constr_lons = iris.Constraint(longitude = lambda cell: West <= cell <= East )
    constr_lats = iris.Constraint(latitude = lambda cell: South <= cell <= North )
    constraints = constr_lats & constr_lons
    
    # stage 3 rough location of perth
    perth3 = slice(400,475,None), slice(400,425,None)
    # stage 5 (smallest nest)
    perth5 = slice(200,1000), slice(200,700) 

    #orog_path = 'data/waroona_oldold/stage3_sfc_orog.nc'
    #orogcubes = fio.read_nc_iris(orog_path)#,constraints=constraints)
    #print(orogcubes)
    ##qplt.pcolormesh(orogcubes[0][0,0], vmax=500)
    #qplt.contourf(orogcubes[0][0,0][perth],40,vmin=-30,vmax=400, cmap='terrain')
    
    
    xwindcubes = fio.read_nc_iris(xwind_path)#,constraints=constraints)
    print("DEBUG: xwindcubes ", xwindcubes)
    xwind = xwindcubes[0][:,:,perth5[0],perth5[1]]
    height = xwind.coord('Hybrid height').points
    #cubetimes = xwind.coord('t')
    #dtimes = utils.dates_from_iris(cubetimes)
    #qplt.contourf(xwind[0,0])
    
    # TODO figure out lats/lons of this run
    latdim = iris.coords.DimCoord(np.linspace(-36,34,xwind.shape[2]),'latitude')
    londim = iris.coords.DimCoord(np.linspace(114,116,xwind.shape[3]),'longitude')
    xwind.add_dim_coord(latdim,2)
    xwind.add_dim_coord(londim,3)
    
    
    
    ywindcubes = fio.read_nc_iris(ywind_path)#,constraints=constraints)
    print("DEBUG: ywindcubes", ywindcubes)
    ywind1 = ywindcubes[0][:,:,perth5[0],perth5[1]]
    ywind1.add_dim_coord(latdim,2) # still staggered
    ywind1.add_dim_coord(londim,3)
    
    zwindcubes = fio.read_nc_iris(zwind_path)#,constraints=constraints)
    print("DEBUG: zwindcubes", zwindcubes)
    zwind = zwindcubes[0][:,:,perth5[0],perth5[1]]
    zwind.add_dim_coord(latdim,2)
    zwind.add_dim_coord(londim,3)

    # font sizes etc
    plotting.init_plots()

    plot_weather_summary(xwind,ywind1,zwind,height,timedim_name='t')

#read_and_plot_oldold_run()
read_and_plot_oldold_run(xwind_path='data/waroona_oldold/combined_alltimes_ml_xwind_stage5.nc',
                        ywind_path='data/waroona_oldold/combined_alltimes_ml_ywind_stage5.nc',
                        zwind_path='data/waroona_oldold/combined_alltimes_ml_zwind_stage5.nc')
