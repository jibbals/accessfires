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

from scipy.interpolate import interp2d

level_constraint = iris.Constraint(model_level_number = lambda cell: cell<=60)
cubes=fio.read_model_run('waroona_run1',datetime(2016,1,6,9), add_z=True,add_winds=True,constraints=level_constraint)

#print (cubes)
start = [-32.84, 115.85]
end = [-32.84, 116.05]
z,w,u, topog = cubes.extract(['z_th','upward_air_velocity','u','surface_altitude'])
lats,lons = w.coord('latitude').points, w.coord('longitude').points
levs = w.coord('model_level_number').points
npoints=50
slicew = utils.cross_section(w[0].data,lats,lons,start,end,npoints=npoints)
slicetopog = utils.cross_section(topog.data,lats,lons,start,end,npoints=npoints)
slicez = utils.cross_section(z.data,lats,lons,start,end,npoints=npoints)
sliceu = utils.cross_section(u[0].data,lats,lons,start,end,npoints=npoints)
xlen=utils.distance_between_points(start,end) # X, Y axes need to be on metres dimension!
xaxis=np.linspace(0,xlen,npoints)
slicex=np.tile(xaxis,(len(levs),1))


x=slicex
y=slicez

plotting.transect_streamplot(x,y,sliceu,slicew)
    