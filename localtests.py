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


# Read waroona
slv,ro1,th1,th2 = fio.read_waroona_iris(dtime0, extent=plotting._extents_['waroona'],add_winds=True)

print(ro1)
s,u,v = ro1.extract(['s','u','v'])
smean = s[:,25:48,:,:].collapsed('model_level_number', iris.analysis.MEAN)
print(smean.level_height)
# metres [z, lat, lon]
zro, = iris.load('data/waroona/umnsaa_2016010515_mdl_ro1.nc', ['height_above_reference_ellipsoid' &
                                                               constr_z & 
                                                               constr_lats & 
                                                               constr_lons])

topog, = fio.read_nc_iris('data/waroona/umnsaa_2016010515_slv.nc',
                         constraints = 'surface_altitude'  & 
                                     constr_lats & 
                                     constr_lons)





