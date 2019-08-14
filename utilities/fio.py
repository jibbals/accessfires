#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:52:09 2019
  
  READING AND WRITING NETCDF AND GRIB(?) FILES
  
@author: jesse
"""

###
## IMPORTS
###

#import numpy as np
from netCDF4 import Dataset
import numpy as np


def read_nc(fpath, keepvars=None):
  '''
    Generic read function for netcdf files
    Reads all dimensions, and all [or some] variables into a dictionary to be returned
  '''
  ncfile =  Dataset(fpath,'r')
  variables = {}
  if keepvars is None:
    for vname in ncfile.variables:
      variables[vname] = ncfile.variables[vname]
  else:
    # keep dimensions and whatever is in keepvars
    for vname in set(keepvars) | set(ncfile.dimensions.keys()):
      variables[vname] = ncfile.variables[vname][:]
  return variables

def read_pcfile(fpath, keepvars=None):
    '''
    Read a umnsaa_pc file, add some extra things such as potential temperature
    '''
    variables=read_nc(fpath,keepvars)
    
    #lat  = variables['latitude']
    #lon  = variables['longitude']
    #lat1 = variables['latitude_0'] # staggered 
    #lon1 = variables['longitude_0'] # staggered
    #z   = variables['model_level_number'  ]
    #z1  = variables['model_level_number_0']
    Ta  = variables['air_temperature'][:,:] # [lat,lon] at surface
    pmsl = variables['air_pressure_at_sea_level'][0,:,:] # [lev,lat,lon]
    #h   = variables['level_height'] # [lev] in metres
    
    # Some variables are on staggered latitude and longitude dimension
    # Some are on staggered time dimension (so there will be 2 time steps instead of none)
    p   = variables['air_pressure'][:,:,:,:] # [time,lev,lat,lon] in pascals
    u1  = variables['x_wind'][:,:,:,:] #[time,levs,lat,lon1] wind speeds are on their directional grid edges
    v1  = variables['y_wind'][:,:,:,:] #[time,levs,lat1,lons]
    #q   = variables['specific_humidity_0'] # [time,lev,lat,lon]
    #w   = variables['upward_air_velocity']# [time,lev,lat,lon] in m/s
    #qc  = variables['mass_fraction_of_cloud_liquid_water_in_air'] + ncfile.variables['mass_fraction_of_cloud_ice_in_air']
    

    # Dimensions
    nt,nz,ny,nx = p.shape

    # height data (z) based on P = P0 exp{ -z/H } with assumed scale height H
    zth = np.zeros(np.shape(p))
    for tstep in range(2):
        zth[tstep] = -(287*300/9.8)*np.log(p[tstep]/pmsl[np.newaxis,:,:])
    
    # Potential temperature based on https://en.wikipedia.org/wiki/Potential_temperature
    # with gas constant R = 287.05 and specific heat capacity c_p = 1004.64
    theta = np.zeros(np.shape(p))
    repTa = np.repeat(Ta[np.newaxis,:,:], nz, axis=0) # repeat Ta along z axis
    for tstep in range(2):
        theta[tstep] = repTa*(1e5/p[tstep])**(287.05/1004.64)

    ## Destagger winds
    #u = np.tile(np.nan,(nz,ny,nx))
    u = np.tile(np.nan,(nt,nz,ny,nx)) # tile repeats the nan accross nz,ny,nx dimensions
    u[:,:,:,1:] = 0.5*(u1[:,:,:,1:] + u1[:,:,:,:-1]) # interpolation of edges
    v = 0.5*(v1[:,:,1::,] + v1[:,:,:-1,:]) # interpolation of edges
    s = np.hypot(u,v) # Speed is hypotenuse of u and v
    # ONE EDGE OF S is now NANs, just set it to adjacent edge speed... (not interested in edge of domain anyway)
    
    variables['theta'] = theta
    variables['z_theta']= zth
    variables['x_wind_destaggered'] = u
    variables['y_wind_destaggered'] = v
    variables['wind_speed'] = s

    
    return variables

def read_topog(pafile='data/umnsaa_pa2016010515.nc'):
    '''
    Read topography and lat,lon from pa file
    '''
    with Dataset(pafile,'r') as ncfile:
        topog = ncfile.variables['surface_altitude'][:,:]
        latt = ncfile.variables['latitude' ][:]
        lont = ncfile.variables['longitude'][:]
    return topog,latt,lont