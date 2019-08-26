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
import timeit

import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt
import matplotlib.patches as mpatches

from utilities import utils,fio,plotting

from finished_plots import clouds_2panel, horizontal_slices

# Try reading file using iris
import iris
import iris.quickplot as qplt


#### SETUP PLOTTING:
plotting.init_plots()

### DATETIMES FOR UM AND FIRE OUTPUT
um_dtimes = np.array([datetime(2016,1,5,15,10) + timedelta(hours=x) for x in range(24)])

dtime0 = um_dtimes[0]
dtimeE = um_dtimes[-1]
um_hour=datetime(dtime0.year,dtime0.month,dtime0.day,dtime0.hour)
ff_dtimes = np.array([um_hour + timedelta(hours=x/60.) for x in range(10,61,10)])


##### TEST CODE HERE


# Constraints on dimensions (testing on local machine, need to save ram)
West,East,South,North = plotting._extents_['waroona']
constr_z = iris.Constraint(model_level_number=lambda cell: cell < 60)
constr_lons = iris.Constraint(longitude = lambda cell: West <= cell <= East )
constr_lats = iris.Constraint(latitude = lambda cell: South <= cell <= North )

# metres [z, lat, lon]
zro, = iris.load('data/waroona/umnsaa_2016010515_mdl_ro1.nc', ['height_above_reference_ellipsoid' &
                                                               constr_z & 
                                                               constr_lats & 
                                                               constr_lons])

topog, = fio.read_nc_iris('data/waroona/umnsaa_2016010515_slv.nc',
                         constraints = 'surface_altitude'  & 
                                     constr_lats & 
                                     constr_lons)

water,ice = fio.read_nc_iris('data/waroona/umnsaa_2016010515_mdl_th2.nc', 
                             constraints = constr_lats & constr_lons).extract(['mass_fraction_of_cloud_liquid_water_in_air',
                                            'mass_fraction_of_cloud_ice_in_air'])



dpi=400
dtime=datetime(2016,1,6,8) # utc 8-10 should be pyrocb
dtime=datetime(2016,1,5,15) # this is output I have on windows
extentname='waroona'
ext='.png'


# Constraints on dimensions (testing on local machine, need to save ram)
West,East,South,North = plotting._extents_[extentname]
constr_z = iris.Constraint(model_level_number=lambda cell: cell < 100)
constr_lons = iris.Constraint(longitude = lambda cell: West <= cell <= East )
constr_lats = iris.Constraint(latitude = lambda cell: South <= cell <= North )
# Read vert motion
_,_,th1cubes,_ = fio.read_waroona_iris(dtime,constraints=constr_z&constr_lats&constr_lons)
w, = th1cubes.extract('upward_air_velocity')
lh = w.coord('level_height').points
lat=w.coord('latitude').points
lon=w.coord('longitude').points

for i in range(6):
    subtime = dtime+timedelta(seconds=600*(i+1))
    horizontal_slices(w[i].data,lh,lat,lon,dtime=subtime)
