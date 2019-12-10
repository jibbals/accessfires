#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker, colors, patches, dates
import numpy as np
import iris

from datetime import datetime, timedelta

from utilities import fio, plotting, utils

_sn_ = 'localtests'

import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from owslib.wmts import WebMapTileService

import cartopy.crs as ccrs
# waroona latlon = -32.84, 115.93
cubes = fio.read_model_run('waroona_run1',fdtime=None,
                           extent=[-32.83,115.92,-32.85,115.94], # just want tiny area at bottom of escarp
                           add_winds=True, add_z=True, add_topog=True,
                           add_RH=True)

x,y,wdir = cubes.extract(['u','v','wind_direction'])
z, topog = cubes.extract(['z_th','surface_altitude'])
T, RH = cubes.extract(['air_temperature','relative_humidity'])

# extract single lat lon
cubes1 = iris.cube.CubeList([x,y,wdir,z,topog, T, RH])
lat,lon = -32.84,115.93
for i in range(len(cubes1)):
    cubes1[i] = cubes1[i].interpolate([('longitude',lon), ('latitude',lat)],
                                      iris.analysis.Linear())

x1,y1,wdir1,z1,topog1,T1,RH1 = cubes1
Tg = T1[:,0].data
RHg = RH1[:,0].data
wdg = wdir1[:,0].data
wsg = np.sqrt(x1[:,0].data**2+y1[:,0].data**2)

# height above ground level
agl = z1.data - topog1.data
agl_coord = iris.coords.AuxCoord(agl, var_name='height_agl', units='m')
agl_coord.guess_bounds()

# interpolate to our heights
heights = np.arange(500,2001,500)
cubes2 = iris.cube.CubeList([x1,y1,wdir1])
## interpolate to 500,1000,1500,2000m above ground level
for i in range(len(cubes2)):
    # add agl coord
    cubes2[i].add_aux_coord(agl_coord,1) 
    cubes2[i] = cubes2[i].interpolate([('height_agl',heights)], iris.analysis.Linear())

print(cubes2)

### Create csv with gridded horizontal wind speed and direction
## 
## Sample headers:
headers=['Time',
         'Temperature',
         'Relative humidity',
         'Wind direction ground',
         'Wind direction 500m',
         'Wind direction 1000m',
         'Wind direction 1500m',
         'Wind direction 2000m',
         'Wind speed ground',
         'Wind speed 500m',
         'Wind speed 1000m',
         'Wind speed 1500m',
         'Wind speed 2000m']
# time format (UTC)
times = utils.dates_from_iris(cubes[0],remove_seconds=True)
tstrings = [dt.strftime("%Y-%m-%dT%H:%M:%S")for dt in times]

import csv
[u2,v2,wd2c] = cubes2
ws2 = np.sqrt(u2.data**2+v2.data**2)
wd2 = wd2c.data
with open('outputs/waroona_winds.csv',mode='w',newline='') as winds_file:
    fwriter=csv.writer(winds_file, delimiter=',')
    fwriter.writerow(headers)
    for i in range(len(times)):
        fwriter.writerow([tstrings[i],Tg[i]-273.15,RHg[i], 
                          wdg[i], wd2[i,0], wd2[i,1], wd2[i,2], wd2[i,3],
                          wsg[i], ws2[i,0], ws2[i,1], ws2[i,2], ws2[i,3]])
fort=None
del fort
