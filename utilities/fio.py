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
from glob import glob
###
## GLOBALS
###


## Sir ivan pc fire files
_topog_sirivan_ = 'data/sirivan/umnsaa_pa2017021121.nc'
_files_sirivan_ = sorted(glob('data/sirivan/umnsaa_pc*.nc'))

def read_nc(fpath, keepvars=None):
  '''
    Generic read function for netcdf files
    Reads all dimensions, and all [or some] variables into a dictionary to be returned
  '''
  ncfile =  Dataset(fpath,'r')
  print("Reading ",fpath, " ... ")
  variables = {}
  if keepvars is None:
    for vname in ncfile.variables:
      variables[vname] = ncfile.variables[vname][:]
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
    variables['zth']= zth
    variables['x_wind_destaggered'] = u
    variables['y_wind_destaggered'] = v
    variables['wind_speed'] = s

    
    return variables

def read_sirivan(fpaths, toplev=80, keepvars=None):
    '''
    Read a subset of the Sir Ivan model outputs
    This will loop over the desired files, so RAM doesn't need to be more than 8GB hopefully
    
    INPUTS:
        fpaths is a list of files to read
        toplev = highest model level to read (save some ram)
        keepvars
            set this if you just one one or two outputs over time
    '''
    # dictionary of data to keep
    data = {}
    
    ## first grab topography from pa vile
    topog,latt,lont = read_topog(_topog_sirivan_)
    data['topog']=topog
    
    ## read first file
    data0 = read_pcfile(fpaths[0])
    
    # Keep dimensions:
    dims = ['latitude','longitude',
            'time_1', # how many time steps in file 
            'model_level_number', # 140 levels
            'pseudo_level', # 5 levels for soil?
            'longitude_0','latitude_0', # staggered lat,lon
            'time_bnds', # file time bounds
            ]
    # copy dims
    for dim in dims:
        data[dim] = data0[dim]
    assert np.all(data['latitude'] == latt), "topog latitude doesn't match pc file latitude"
    assert np.all(data['longitude'] == lont), "topog longitude doesn't match pc file longitude"
    
    ## data without time dim
    flatdata = ['air_temperature',      # [lat, lon] K at surface
                'specific_humidity',    # [lat,lon] kg/kg
                'surface_air_pressure', # [lat,lon] Pa
                'surface_temperature',  # [lat,lon] K
                #'level_height',         # [z] metres from ground level]
                ]
    
    ## 4 dim data [t1, z, lat, lon]
    fulldata = ['air_pressure',                                 # Pascals
                'air_temperature_0',                            # K
                'mass_fraction_of_cloud_ice_in_air',            # kg/kg
                'mass_fraction_of_cloud_liquid_water_in_air',   # kg/kg
                'theta',                                        # K (potential temp)
                'wind_speed',                                   # m/s (wind speed)
                'specific_humidity_0',                          # kg/kg
                'upward_air_velocity',                          # m/s
                'x_wind_destaggered',                           # m/s
                'y_wind_destaggered',                           # m/s
                'zth',                                          # m (height approximated)
                ]
    if keepvars is not None:
        fulldata = list(set(fulldata) & set(keepvars))
    
    ## Start tracking the time dimensions
    hour = []
    time = []
    
    # Gregorian hours since 19700101 00:00:00
    hour.append(data0['time_0']) 
    time.append(data0['time_1'][0])
    time.append(data0['time_1'][1])
    
    ## Copy all the data we want
    for varname in flatdata:
        # add time dim to the flat data
        data[varname] = data0[varname][np.newaxis]
    
    for varname in fulldata:
		# Just read up to toplev 
        data[varname] = data0[varname][:,:toplev,:,:]
        
    ## also copy height
    data['level_height'] = data0['level_height'][:toplev]
    
    del data0
    
    # read each file and append along time dimension
    for i,fpath in enumerate(fpaths):
        if i==0:
            continue
        datai = read_pcfile(fpath)
        
        # Append each file along the time dimension
        for varname in flatdata+fulldata:
            data[varname] = np.append(data[varname], datai[varname], axis=0)
        
        # Append time step
        hour.append(datai['time_0']) 
        time.append(datai['time_1'][0])
        time.append(datai['time_1'][1])
        
        del datai
    
    data['hour'] = np.array(hour)
    data['time'] = np.array(time)
    # also for convenience add cloud metric
    if keepvars is None:
        data['qc'] = data['mass_fraction_of_cloud_ice_in_air'] + data['mass_fraction_of_cloud_liquid_water_in_air']
    return data

def read_topog(pafile='data/umnsaa_pa2016010515.nc'):
    '''
    Read topography and lat,lon from pa file
    '''
    with Dataset(pafile,'r') as ncfile:
        topog = ncfile.variables['surface_altitude'][:,:]
        latt = ncfile.variables['latitude' ][:]
        lont = ncfile.variables['longitude'][:]
    return topog,latt,lont
