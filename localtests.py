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

from finished_plots import clouds_2panel

# Try reading file using iris
import iris
import iris.quickplot as qplt
from iris.fileformats.netcdf import load_cubes

plotting.init_plots()

um_dtimes = np.array([datetime(2016,1,5,15,10) + timedelta(hours=x) for x in range(24)])

dtime = um_dtimes[0]
extentname='waroona'
vectorskip=9


# Constraints on dimensions
West,East,South,North = plotting._extents_['waroona']
constr_z = iris.Constraint(model_level_number=lambda cell: cell < 60)
constr_lons = iris.Constraint(longitude = lambda cell: West <= cell <= East )
constr_lats = iris.Constraint(latitude = lambda cell: South <= cell <= North )

# metres [z, lat, lon]
zro, = iris.load('data/waroona/umnsaa_2016010515_mdl_ro1.nc', ['height_above_reference_ellipsoid' &
                                                               constr_z & 
                                                               constr_lats & 
                                                               constr_lons])

topography = fio.read_nc_iris('data/waroona/umnsaa_2016010515_slv.nc',
                              constraints = 'surface_altitude'  & 
                                            constr_lats & 
                                            constr_lons)[0]


# Read the cubes
slv,ro1,th1,th2 = fio.read_waroona_iris(dtime=dtime, 
                                        constraints = [constr_z &
                                                       constr_lons &
                                                       constr_lats])

# wind speeds need to be interpolated onto non-staggered latlons
# u is on lon1 dim
# v is on lat1 dim

# pull out wind speed
p, u1, v1 = ro1.extract(['air_pressure','x_wind','y_wind'])
# DESTAGGER u and v using iris interpolate
# u1: [t,z, lat, lon1]
# v1: [t,z, lat1, lon]  # put these both onto [t,z,lat,lon]
u = u1.interpolate([('longitude',p.coord('longitude').points)],
                   iris.analysis.Linear())
v = v1.interpolate([('latitude',p.coord('latitude').points)],
                   iris.analysis.Linear())
lon=u.coord('longitude').points
lat=u.coord('latitude').points
lon1=u1.coord('longitude').points
lat1=v1.coord('latitude').points

# Get wind speed cube
s = iris.analysis.maths.apply_ufunc(np.hypot,u,v)

sh,Ta = th1.extract(['specific_humidity','air_temperature'])

rh = utils.relative_humidity_from_specific(sh.data,Ta.data)

qc1,qc2 = th2.extract(['mass_fraction_of_cloud_ice_in_air','mass_fraction_of_cloud_liquid_water_in_air'])

qc = qc1+qc2


# datetime of outputs
tdim = p.coord('time')

d0 = datetime.strptime(str(p.units),'seconds since %Y-%m-%d %H:%M:00')
timesteps = utils.date_from_gregorian(tdim.points, d0=do)

# also loop over different transects
for i_transect in np.arange(1,6.5,1, dtype=int):
    for tstep in range(len(timesteps)):
        clouds_2panel(
                      dtime=timesteps[tstep],
                      extentname=extentname,
                      vectorskip=vectorskip, 
                      transect=i_transect)




