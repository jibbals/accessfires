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

#import sys

from utilities import utils, plotting, fio

plotting.init_plots()
extents = plotting._extents_
latlons = plotting._latlons_

# Data file locations
pafile = 'data/umnsaa_pa2016010515.nc'
pcfile = 'data/umnsaa_pc2016010515.nc'
with Dataset(pafile,'r') as ncfile:
    topog = ncfile.variables['surface_altitude'][:,:]
    latt = ncfile.variables['latitude' ][:]
    lont = ncfile.variables['longitude'][:]

plt.close()


# Google map image tiles view of synoptic map
fig,ax,proj=plotting.map_google(extents['waroonas'],
                                zoom=6,
                                subplotxyn=[2,1,1],
                                gridlines=[np.arange(-50,-10,2),
                                           np.arange(100,150,4)])

## Add box around zoomed in area
ax.add_patch(mpatches.Rectangle(xy=latlons['waroona'][::-1], 
                                width=.4, 
                                height=.3,
                                facecolor='red', 
                                alpha=0.4, 
                                transform=ccrs.PlateCarree()))
## add text?


## Add scale
scaleloc=(0.2,0.05)
plotting.scale_bar(ax,proj,100, location=scaleloc)


## Look at waroona and yarloop
_,ax2,_ = plotting.map_google(extents['waroona'],zoom=10,fig=fig,subplotxyn=[2,2,3],draw_gridlines=False)

## Add scale
plotting.scale_bar(ax2,proj,10, location=scaleloc)

## Add contour plot showing topography
plt.subplot(2,2,4)
plt.contourf(lont,latt,topog, cmap='copper')
# set x and y limits to match extent
xlims = extents['waroona'][0:2] # East to West
ylims = extents['waroona'][2:] # South to North
plt.ylim(ylims); plt.xlim(xlims)
plt.colorbar(label="m")
plt.title("Topography")
## Turn off the tick values
plt.xticks([]); plt.yticks([])
plt.show()

assert False,"stop here"




nlev=60 # Save on ram just look at lower levels
ztop=3000
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

# wind speed
s[np.isnan(s)] = -5000 # There is one row or column of s that is np.NaN, one of the edges I think
slices = utils.cross_section(s,start,end,lat,lon,npoints=npoints)
slevels = np.arange(0,22,2)
cmaps = plt.cm.get_cmap('YlGnBu')

# clouds by water+ice content
sliceqc = utils.cross_section(qc,start,end,lat,lon,npoints=npoints)
qccontours = np.arange(0.0,2.25,0.25)
cmapqc = plt.cm.get_cmap('YlGnBu')

for ax, slicedata, slicelevels, cmap, slicecontours, title,norm, cbarform in [
    [axes[0,0], slicetheta,thetalevels,thetacmap,thetacontours, "T$_{\\theta}$ (K)", None, None],
    [axes[0,1], slicew, wlevels, cmapw, wcontours, "Vertical motion (m/s)",wnorm, tick.ScalarFormatter()],
    [axes[1,0], slices, slevels, cmaps, slevels, "Wind (m/s)", None, None],
    [axes[1,1], sliceqc, qccontours, cmapqc, qccontours, "Water+ice (kg/kg air)",None,None]
    ]:
  plt.sca(ax)
  # Note that contourf can work with non-plaid coordinate grids provided both are 2-d
  # Contour inputs: xaxis, yaxis, data, colour gradient 
  plt.contourf(slicex,slicez,slicedata,slicelevels,cmap=cmap,norm=norm)
  plt.colorbar(format=cbarform)
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
  plt.xlabel("Transect")
  #plt.xticks(xticks,xlabels, rotation=25)
  plt.xticks(None,None)
# Pull out cross section of potential temperature

