# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 10:50:05 2019
    Example for reading fire output
@author: Mika Peace
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import time
from datetime import datetime, timedelta
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy import interpolate
from numpy import ma

with PdfPages('multipage_pdf.pdf') as pdf:

    datadir = '/short/en0/hxy548/tmp/waroona/0p3/'

    txr=range(0,780,5)
   # txr=range(1000,1440,5)
    for tx in txr:
        t0=0

        ncfile = Dataset(datadir+'f1_terrain_height.CSIRO_24h.nc', 'r')
        lon = ncfile.variables['lon'][:]
        lat = ncfile.variables['lat'][:]
        terrain = ncfile.variables['terrain_height'][t0,:,:].transpose()
        ncfile.close()

        ncfile = Dataset(datadir+'sensible_heat.CSIRO_24h.20160107T0300Z.nc', 'r')
        sheat = ncfile.variables['SHEAT_2'][tx,:,:].transpose()
        time = ncfile.variables['time'][:]
        ncfile.close()

        ncfile = Dataset(datadir+'fire_speed.CSIRO_24h.20160107T0300Z.nc', 'r')
        firespeed = ncfile.variables['fire_speed'][tx,:,:].transpose()
        time = ncfile.variables['time'][:]
        ncfile.close()

        ncfile = Dataset(datadir+'firefront.CSIRO_24h.20160107T0300Z.nc', 'r')
        ff = ncfile.variables['firefront'][tx,:,:].transpose()
        time = ncfile.variables['time'][tx]
        ncfile.close()

        ncfile = Dataset(datadir+'10m_uwind.CSIRO_24h.20160107T0300Z.nc', 'r')
        u = ncfile.variables['UWIND_2'][tx,:,:].transpose()
        u1 = u[::8,::8]
        ncfile.close()

        ncfile = Dataset(datadir+'10m_vwind.CSIRO_24h.20160107T0300Z.nc', 'r')
        v = ncfile.variables['VWIND_2'][tx,:,:].transpose()
        v1 = v[::8,::8]
        ncfile.close()

        wspd=np.hypot(u,v)#wind speed
        wdir=(270-(np.arctan2(v,u)*180/np.pi))%360
        ts=datetime(2016,01,07,03,00,00)#starttime from netcdf file
        mytime=(ts+timedelta(seconds=time))#time adds to deconds from .nc file
        print(mytime)

        lon1 = lon[::8]
        lat1 = lat[::8]
        lons,lats = np.meshgrid(lon,lat)
        lons1,lats1 = np.meshgrid(lon1,lat1)
        #fig = plt.figure()
        plt.figure()
        plt.clf()

        plt.subplot(221,aspect=1.0)
        clevs= np.linspace(-150,550,29,endpoint=True)#linspace seems to work better than arange in colorbar
        cmaptr=plt.cm.get_cmap("terrain")
        plt.contourf(lons,lats,terrain,clevs,cmap=cmaptr)
        cb=plt.colorbar(ticks=clevs, fraction=0.045, pad=0.01)
        cb.set_label(' ', size=20)
        cb.ax.tick_params(labelsize=6)
        plt.contour(lons,lats,ff,np.array([0.0]), colors='red')

        ##annotations
        plt.annotate('Fire ignition', xy=(116.2,-32.9), xytext=(116.10,-32.92), fontsize=6)
        plt.plot([116.17],[-32.89], 'ro', ms=3)
        plt.annotate('Waroona', xy=(116.2,-32.8), xytext=(115.9,-32.82), fontsize=6)
        plt.plot([115.93],[-32.84], 'ko', ms=3)
        plt.annotate('Yarloop', xy=(116.1,-33.1), xytext=(115.85,-32.99), fontsize=6)
        plt.plot([115.90],[-32.96], 'ko', ms=3)
        ##sets limits and increments on x and y axes
        plt.xlim([115.6,116.2])# was [115.8,116.2]
        plt.ylim([-33.1,-32.7])# was [-33.1,32.7]
        lonlab = (np.arange(115.6,116.2,0.1)) #was [115.8,116.2,0.1]
        latlab = (np.arange(-33.1,-32.7,0.1)) #was [33.1-32.7,0.1]
        plt.xticks(lonlab,fontsize=6)
        plt.yticks(latlab,fontsize=6)
        plt.ticklabel_format(useOffset=False)#keeps the tick lables from switching to scientific format
        plt.title('Topography and fire perimeter', fontsize=8)

        plt.subplot(222,aspect=1.0)
        clevs= np.linspace(0,5.5,12, endpoint=True)#linspace seems to work beter than arange in colorbar
        cmap2=plt.cm.get_cmap("Oranges")
        plt.contourf(lons,lats,np.log10(sheat+1),clevs,cmap=cmap2)
        cb=plt.colorbar(ticks=clevs, fraction=0.045, pad=0.01)
        cb.set_label(' ', size=20)
        cb.ax.tick_params(labelsize=6)
        cmap2=plt.cm.get_cmap("YlOrRd")

        ##annotations
        plt.annotate('Fire ignition', xy=(116.2,-32.9), xytext=(116.10,-32.92), fontsize=6)
        plt.plot([116.17],[-32.89], 'ro', ms=3)
        plt.annotate('Waroona', xy=(116.2,-32.8), xytext=(115.9,-32.82), fontsize=6)
        plt.plot([115.93],[-32.84], 'ko', ms=3)
        plt.annotate('Yarloop', xy=(116.1,-33.1), xytext=(115.85,-32.99), fontsize=6)
        plt.plot([115.90],[-32.96], 'ko', ms=3)
        ##sets limits and increments on x and y axes
        plt.xlim([115.6,116.2])# was [115.8,116.2]
        plt.ylim([-33.1,-32.7])# was [-33.1,32.7]
        lonlab = (np.arange(115.6,116.2,0.1)) #was [115.8,116.2,0.1]
        latlab = (np.arange(-33.1,-32.7,0.1)) #was [33.1-32.7,0.1]
        plt.xticks(lonlab,fontsize=6)
        plt.yticks(latlab,fontsize=6)
        plt.ticklabel_format(useOffset=False)#keeps the tick lables from switching to scientific format
        plt.title('log(Sens heat+1)', fontsize=8)

        plt.subplot(223,aspect=1.0)
        #contour levels for colorbar
        clevs= np.linspace(0,20,11, endpoint=True)#linspace seems to work beter than arange in colorbar
        cmap1=plt.cm.get_cmap("PuBuGn") # was "summer"
        plt.contourf(lons,lats,wspd,clevs,cmap=cmap1)
        cb=plt.colorbar(ticks=clevs, fraction=0.045, pad=0.01)
        cb.set_label(' ', size=20)
        cb.ax.tick_params(labelsize=6)

        ##annotations
        plt.annotate('Fire ignition', xy=(116.2,-32.9), xytext=(116.10,-32.92), fontsize=6)
        plt.plot([116.17],[-32.89], 'ro', ms=3)
        plt.annotate('Waroona', xy=(116.2,-32.8), xytext=(115.9,-32.82), fontsize=6)
        plt.plot([115.93],[-32.84], 'ko', ms=3)
        plt.annotate('Yarloop', xy=(116.1,-33.1), xytext=(115.85,-32.99), fontsize=6)
        plt.plot([115.90],[-32.96], 'ko', ms=3)
        ##sets limits and increments on x and y axes
        plt.xlim([115.6,116.2])# was [115.8,116.2]
        plt.ylim([-33.1,-32.7])# was [-33.1,32.7]
        lonlab = (np.arange(115.6,116.2,0.1)) #was [115.8,116.2,0.1]
        latlab = (np.arange(-33.1,-32.7,0.1)) #was [33.1-32.7,0.1]
        plt.xticks(lonlab,fontsize=6)
        plt.yticks(latlab,fontsize=6)
        plt.ticklabel_format(useOffset=False)#keeps the tick lables from switching to scientific format
        plt.title('Wind speed (m/s)', fontsize=8)
    #
        plt.subplot(224,aspect=1.0)
        clevs= np.linspace(0,360,25, endpoint=True)#linspace seems to work beter than arange in colorbar
        cmap2=plt.cm.get_cmap("rainbow")
        plt.contourf(lons,lats,wdir,clevs,cmap=cmap2)
        cb=plt.colorbar(ticks=clevs, fraction=0.045, pad=0.01)
        cb.set_label(' ', size=20)
        cb.ax.tick_params(labelsize=6)
        plt.contour(lons,lats,ff,np.array([0.0]), colors='red')
    #
        ##annotations
        plt.annotate('Fire ignition', xy=(116.2,-32.9), xytext=(116.10,-32.92), fontsize=6)
        plt.plot([116.17],[-32.89], 'ro', ms=3)
        plt.annotate('Waroona', xy=(116.2,-32.8), xytext=(115.9,-32.82), fontsize=6)
        plt.plot([115.93],[-32.84], 'ko', ms=3)
        plt.annotate('Yarloop', xy=(116.1,-33.1), xytext=(115.85,-32.99), fontsize=6)
        plt.plot([115.90],[-32.96], 'ko', ms=3)
        ##sets limits and increments on x and y axes
        plt.xlim([115.6,116.2])# was [115.8,116.2]
        plt.ylim([-33.1,-32.7])# was [-33.1,32.7]
        lonlab = (np.arange(115.6,116.2,0.1)) #was [115.8,116.2,0.1]
        latlab = (np.arange(-33.1,-32.7,0.1)) #was [33.1-32.7,0.1]
        plt.xticks(lonlab,fontsize=6)
        plt.yticks(latlab,fontsize=6)
        plt.ticklabel_format(useOffset=False)#keeps the tick lables from switching to scientific format
        plt.title('Wind direction (degN)', fontsize=8)
    #    #show plot
        plt.suptitle(str(mytime))
    #    plt.savefig(str(mytime)+'wind.png')
        plt.savefig(str(mytime)+'wind.png', dpi=600)
        plt.close()
        plt.cla()

