#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""

import cartopy
import cartopy.io.img_tiles as cimgt
import matplotlib as mpl
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
from netCDF4 import Dataset
#import sys

from utilities import utils, plotting, fio

plotting.init_plots()
extents = plotting._extents_
latlons = plotting._latlons_

plt.close()

# Data file locations
pafile = 'data/umnsaa_pa2016010515.nc'
pcfile = 'data/umnsaa_pc2016010515.nc'

nlev=60 # Save on ram just look at lower levels
ztop=4000
with Dataset(pcfile,'r') as ncfile:
    
    ## PULL OUT MASKED ARRAYS
    #zth  = ncfile.variables['height_theta'][0,:,::-1,:] # flip latitudes for interpolate
    #zrho = ncfile.variables['height_rho'  ][0,:,::-1,:]
    lat  = ncfile.variables['latitude'   ][:]
    lon  = ncfile.variables['longitude'  ][:]
    lat1 = ncfile.variables['latitude_0' ][:] # Edges??
    lon1 = ncfile.variables['longitude_0'][:] # Edges??
    z   = ncfile.variables['model_level_number'  ][:]
    z1  = ncfile.variables['model_level_number_0'][:]
    Ta  = ncfile.variables['air_temperature'][:,:] # at surface? only stored at surface?
    p   = ncfile.variables['air_pressure'][0,0:nlev,:,:] # in pascals
    pmsl = ncfile.variables['air_pressure_at_sea_level'][0,:,:]
    u1  = ncfile.variables['x_wind'][0,0:nlev,:,:] # wind speeds are on their directional grid edges (z,lats,lon1)
    v1  = ncfile.variables['y_wind'][0,0:nlev,:,:] # [z,lat1,lons]
    q   = ncfile.variables['specific_humidity_0'][0,0:nlev,:,:]
    h   = ncfile.variables['level_height'][:] # in metres
    w   = ncfile.variables['upward_air_velocity'][0,0:nlev,:,:]
    qc  = ncfile.variables['mass_fraction_of_cloud_liquid_water_in_air'][0,0:nlev,:,:] + ncfile.variables['mass_fraction_of_cloud_ice_in_air'][0,0:nlev,:,:]
    

# Dimensions
nz,ny,nx = p.shape

# Fudge some height data
zth = -(287*300/9.8)*np.log(p/pmsl[np.newaxis,:,:])
zrho = zth
theta = Ta*(1e5/p)**(287.05/1004.64)

## Destagger winds
#u = np.tile(np.nan,(nz,ny,nx))
u = np.tile(np.nan,(nz,ny,nx)) # tile repeats the nan accross nz,ny,nx dimensions
u[:,:,1:] = 0.5*(u1[:,:,1:] + u1[:,:,:-1]) # interpolation of edges
v = 0.5*(v1[:,1::,] + v1[:,:-1,:]) # interpolation of edges
s = np.hypot(u,v) # Speed is hypotenuse of u and v
lonu = lon1
latu = lat
lonv = lon
latv = lat1
    


# dummy topog for now
topog = np.zeros([len(lat),len(lon)])
latt,lont = lat,lon
# READ TOPOG DATA FROM PA
with Dataset(pafile,'r') as ncfile:
    topog = ncfile.variables['surface_altitude'][:,:]
    latt = ncfile.variables['latitude' ][:]
    lont = ncfile.variables['longitude'][:]



# cross section lat,lon start and finish
start,end = [-32.75,115.75],[-32.95,116]
npoints = 40

# Pull out cross section of topography and height
slicetopog = utils.cross_section(topog,start,end, latt,lont,npoints=npoints)
slicez = utils.cross_section(zth,start,end,lat,lon,npoints=npoints)
xticks,xlabels = utils.cross_section_ticks_labels(start,end)
xaxis=np.linspace(0,1,npoints)
slicex=np.tile(xaxis,(nz,1))

plt.figure()
plt.plot(xaxis,slicetopog)
plt.xticks(xticks,xlabels)
plt.ylabel('Height (m)')
#plt.show()


## First image shows potential temp vs height
f,axes = plt.subplots(2,2, sharex=True, sharey=True)

# Potential temperature
slicetheta = utils.cross_section(theta,start,end,lat,lon,npoints=npoints)
thetalevels = np.arange(280,320,2) # contour lines to plot
thetacontours = thetalevels
thetacmap = plt.cm.get_cmap('YlOrRd')

# Vertical velocity
slicew = utils.cross_section(w,start,end,lat,lon,npoints=npoints)
wlevels = 2.0**np.arange(-2,6)
wlevels = np.union1d(np.union1d(wlevels,-wlevels),np.array([0]))
wcontours = np.array([0])
cmapw = plt.cm.get_cmap('PiYG')
cmapw.set_over('k')
wnorm = col.SymLogNorm(0.25)

for ax, slicedata, slicelevels, cmap, slicecontours, title,norm in [
    [axes[0,0], slicetheta,thetalevels,thetacmap,thetacontours, "T$_{\theta}$ (K)", None],
    axes[0,1], slicew, wlevels, cmapw, wcontours, "Vertical motion (m/s)",wnorm]:
  plt.sca(ax)
  # Note that contourf can work with non-plaid coordinate grids provided both are 2-d
  # Contour inputs: xaxis, yaxis, data, colour gradient 
  plt.contourf(slicex,slicez,slicedata,slicelevels,cmap=cmap,norm=norm)
  plt.colorbar()
  # Add contour lines
  plt.contour(slicex,slicez,slicedata,slicecontours,colors='k')            
  # make sure land is obvious
  plt.fill_between(xaxis,slicetopog,interpolate=True,facecolor='black')
  plt.xticks(None,None)
  if ztop != None:
    plt.ylim(0,ztop)
  plt.xlabel('')
  plt.title(title)

for ax in [axes[0,0],axes[1,0]]:
  plt.sca(ax)
  plt.ylabel('Height (m)')

for ax in [axes[1,0],axes[1,1]]:
  plt.sca(ax)
  #plt.xticks(xticks,xlabels,rotation=-15)
  plt.xlabel("Horizontal")
# Pull out cross section of potential temperature

