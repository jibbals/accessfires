#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""


# IMPORTS

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import iris # so we can add level number constraint
import iris.analysis
import iris.quickplot as qplt

# read tiff file
from osgeo import gdal, osr
import cartopy.crs as ccrs

from utilities import utils, plotting, fio, constants

from scipy import interpolate

#print (cubes)
start = [-32.84, 115.85]
end = [-32.84, 116.05]
fdtime = datetime(2016,1,6,12)
mr='waroona_run1'
import pyrocb

#pyrocb.pyrocb_model_run(model_run=mr,dtime=fdtime)
#pyrocb.moving_pyrocb(model_run=mr,dtimes=[fdtime,])


level_constraint = iris.Constraint(model_level_number = lambda cell: cell<=60)
cubes=fio.read_model_run('waroona_run1',datetime(2016,1,6,9), add_z=True,add_winds=True,constraints=level_constraint)
z, w, u, topog = cubes.extract(['z_th','upward_air_velocity','u','surface_altitude'])
lats,lons = w.coord('latitude').points, w.coord('longitude').points
nz = z.shape[0]
nx = utils.number_of_interp_points(lats,lons,start,end)
lataxis,lonaxis = utils.latslons_axes_along_transect(lats,lons,start,end,nx)
z_th = z.data
# Interpolate data along slice (in model height coordinate). Note that the
# transpose is needed because RectBivariateSpline assumes the axis order (x,y)
slicedata = np.tile(np.nan, [nz, nx])
slicez = np.tile(np.nan, [nz, nx])
data=w[0].data
for k in range(0,nz):
    f = interpolate.RectBivariateSpline(lons,lats,data[k,:,:].T)
    slicedata[k,:] = f.ev(lonaxis,lataxis)
    if z_th is not None:
        f2 = interpolate.RectBivariateSpline(lons,lats,z_th[k,:,:].T)
        slicez[k,:] = f2.ev(lonaxis, lataxis)

#cubes=fio.read_model_run('waroona_run1',datetime(2016,1,6,9), add_z=True,add_winds=True,constraints=level_constraint)
#z,w,u, topog = cubes.extract(['z_th','upward_air_velocity','u','surface_altitude'])
#
#levs = w.coord('model_level_number').points
#npoints=50
#slicew = utils.cross_section(w[0].data,lats,lons,start,end,npoints=npoints)
#slicetopog = utils.cross_section(topog.data,lats,lons,start,end,npoints=npoints)
#slicez = utils.cross_section(z.data,lats,lons,start,end,npoints=npoints)
#sliceu = utils.cross_section(u[0].data,lats,lons,start,end,npoints=npoints)
#xlen=utils.distance_between_points(start,end) # X, Y axes need to be on metres dimension!
#xaxis=np.linspace(0,xlen,npoints)
#slicex=np.tile(xaxis,(len(levs),1))
#
#
#x=slicex
#y=slicez
#
#plotting.transect_streamplot(x,y,sliceu,slicew)
    