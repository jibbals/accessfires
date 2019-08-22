#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""

import matplotlib as mpl
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
from netCDF4 import Dataset

from datetime import datetime, timedelta

import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt
import matplotlib.patches as mpatches

from utilities import utils,fio,plotting

# Try reading file using iris
import iris
import iris.quickplot as qplt
from iris.fileformats.netcdf import load_cubes

plotting.init_plots()

fpath='C:/Users/jgreensl/Desktop/data/waroona_fire/firefront.CSIRO_24h.20160105T1500Z.nc'

#ff, = fio.read_nc_iris(fpath)
#ff = ff[9::10] # take every 10th minute

ff, t, lats,lons = fio.read_fire(fpath)
dtimes=utils.date_from_gregorian(t/3600.,d0=datetime(2016,1,5,15))

topog,latt,lont = fio.read_topog()
clevs= np.linspace(-150,550,50,endpoint=True)#linspace seems to work better than arange in colorbar

plt.close()
plt.figure()
for i, j in enumerate(np.arange(1,144+1,144/6.0, dtype=int)):
    plt.subplot(3,2,i+1)
    print(i,j)
    print(ff[j-1].shape, dtimes[j-1])
    cmaptr=plt.cm.get_cmap("terrain")    
    plt.contourf(lons,lats,topog,levels=clevs,cmap=cmaptr)
    plt.contour(lons,lats,np.transpose(ff[j-1]),np.array([0]), colors='red')
    plt.title(dtimes[j-1])


