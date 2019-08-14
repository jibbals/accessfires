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

import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt
import matplotlib.patches as mpatches

from utilities import utils,fio,plotting

dat = fio.read_pcfile('data/umnsaa_pc2016010515.nc')

nlev=80 # Save on ram just look at lower levels
tstep=0

lat  = dat['latitude'   ][:]
lon  = dat['longitude'  ][:]
lat1 = dat['latitude_0' ][:] # Staggered
lon1 = dat['longitude_0'][:] # Staggered
z   = dat['model_level_number'  ][:]
z1  = dat['model_level_number_0'][:]
Ta  = dat['air_temperature'][:,:] # at surface
p   = dat['air_pressure'][tstep,:nlev,:,:] # in pascals
pmsl = dat['air_pressure_at_sea_level'][0,:,:] # just take surface 
u1  = dat['x_wind'][tstep,:nlev,:,:] # wind speeds are on their directional grid edges (z,lats,lon1)
v1  = dat['y_wind'][tstep,:nlev,:,:] # [z,lat1,lons]
q   = dat['specific_humidity_0'][tstep,:nlev,:,:]
lh   = dat['level_height'][:nlev] # in metres (grid box middle, distance from ground)
w   = dat['upward_air_velocity'][tstep,:nlev,:,:]
qc  = dat['mass_fraction_of_cloud_liquid_water_in_air'][tstep,:nlev,:,:] + dat['mass_fraction_of_cloud_ice_in_air'][tstep,:nlev,:,:]
# Also some stuff based on calculated data (in fio.py)
theta = dat['theta'][tstep,:nlev,:,:]
zth = dat['z_theta'][tstep,:nlev,:,:]
u   = dat['x_wind_destaggered'][tstep,:nlev,:,:]
v   = dat['y_wind_destaggered'][tstep,:nlev,:,:]
s   = dat['wind_speed'][tstep,:nlev,:,:]
# Save some ram:
del dat

# Dimensions
nz,ny,nx = p.shape

# READ TOPOG DATA FROM PA
topog, latt, lont = fio.read_topog()
# 400x400 (not matching the pc file)

# make actual height data
replh = np.repeat(lh[:,np.newaxis],400,axis=1)
replh = np.repeat(replh[:,:,np.newaxis],400,axis=2)
h  = topog[np.newaxis,:,:] + replh

# Compare height (m) to zth (m), should be similar
# need to match lat/lon for proper look
plt.scatter(h[:,50,50],zth[:,50,50], )
plt.plot([0,h[-1,50,50]],[0,h[-1,50,50]],'k--',alpha=0.6,linewidth=2,label='1:1')
plt.xlabel("model level height + topog")
plt.ylabel("zth")
plt.legend()
plt.title("h vs zth, not index matched (yet)")

#assert False,"stop here"


# cross section lat,lon start and finish
transects = plotting._transects_
start,end = transects['waroona1']
npoints = 100

# Pull out cross section of topography and height
slicetopog = utils.cross_section(topog,latt,lont,start,end,npoints=npoints)
xticks,xlabels = utils.cross_section_ticks_labels(start,end)
xaxis=np.linspace(0,1,npoints)
slicex=np.tile(xaxis,(nz,1))

plt.figure()
plt.plot(xaxis,slicetopog)
plt.xticks(xticks,xlabels)
plt.ylabel('Height (m)')
#plt.show()


## Image showing topography (with wind speeds?) and transects of vert motion and wind speed
plt.close()
waroona=plotting._extents_['waroona']
f = plt.figure(figsize=[7,10])

ax1 = plt.subplot(2,1,1)
plotting.map_topography(waroona,topog,latt,lont)
plt.title('Topography, wind speed')
# start to end x=[lon0,lon1], y=[lat0, lat1]
plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
         linewidth=2, marker='o', markersize=5, markercolor='red')

ax2 = plt.subplot(2,2,3)
plotting.transect_w(w,zth,lat,lon,start,end,topog=topog,latt=latt,lont=lont)

ax3 = plt.subplot(2,2,4)
plotting.transect_s(s,zth,lat,lon,start,end,topog=topog,latt=latt,lont=lont)
plt.yticks([])

plt.show()