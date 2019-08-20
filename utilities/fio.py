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
_topog_waroona_ = 'data/waroona/umnsaa_pa2016010515.nc'
_files_sirivan_ = sorted(glob('data/sirivan/umnsaa_pc*.nc'))
# Old version of waroona output:
_files_waroona_old_= sorted(glob('data/waroona/umnsaa_pc*.nc'))

# 3 file types for new waroona output
## slx - single level output
## mdl_ro1 - multi level, on ro dimension (winds)
## mdl_th1 - multi level, on theta dimension
## mdl_th2 - multi level, on theta dimension (since um output file size is limited)
_file_types_waroona_ = {'slv':['specific_humidity',  # [t,lat,lon]
                               'surface_air_pressure',
                               'surface_altitude', # [lat,lon] # Is this the new topog? 
                               'surface_temperature', # [t,y,x]
                               'time', # = 6 ;
                               #'level_height', # = 140 ;
                               'latitude', # = 576 ;
                               'longitude', # = 576
                               # x_wind and y_wind are here also...
                               ],
                        'mdl_ro1':['air_pressure', # [t,z,lat,lon]
                               'x_wind', # [t,z,y,x]
                               'y_wind', # [t,z,y,x]
                               'height_above_reference_ellipsoid', #[z,y,x] # new level heights?
                               'time', # = 6 ;
                               'level_height', # = 140 ;
                               'latitude', # = 576 ;
                               'longitude', # = 576
                               ],
                        'mdl_th1':['air_pressure', # Pa [t,z,y,x]
                                   'air_temperature', # K [t,z,y,x]
                                   'height_above_reference_ellipsoid', # m [t,z,y,x]
                                   'specific_humidity', # kg/kg [tzyx]
                                   'upward_air_velocity', # m/s [tzyx]
                                   'time', # = 6 ;
                                   'level_height', # = 140 ;
                                   'latitude', # = 576 ;
                                   'longitude', # = 576
                                   ],
                        'mdl_th2':['mass_fraction_of_cloud_ice_in_air', # kg/kg [tzyx]
                                   'mass_fraction_of_cloud_liquid_water_in_air',
                                   'time', # = 6 ;
                                   'level_height', # = 140 ;
                                   'latitude', # = 576 ;
                                   'longitude', # = 576
                                   ],
                       }
_files_waroona_ = 'data/waroona/umnsaa_%s_%s.nc' # umnsaa_YYYYMMDDhh_slx.nc


def read_nc(fpath, keepvars=None):
  '''
    Generic read function for netcdf files
    Reads all dimensions, and all [or some] variables into a dictionary to be returned
  '''
  ncfile =  Dataset(fpath,'r')
  print("INFO: Reading ",fpath, " ... ")
  variables = {}
  if keepvars is None:
    for vname in ncfile.variables:
      variables[vname] = ncfile.variables[vname][:]
  else:
    # keep dimensions and whatever is in keepvars
    for vname in set(keepvars) | (set(ncfile.dimensions.keys()) & set(ncfile.variables.keys())):
      variables[vname] = ncfile.variables[vname][:]

  print("INFO: finished reading ",fpath)
  return variables

def read_pcfile(fpath, keepvars=None):
    '''
    Read a umnsaa_pc file, add some extra things such as potential temperature
    '''
    variables=read_nc(fpath,keepvars)
    
    if 'air_pressure' in variables.keys():
        p   = variables['air_pressure'][:,:,:,:] # [time,lev,lat,lon] in pascals
        nt,nz,ny,nx = p.shape

    # height data (z) based on P = P0 exp{ -z/H } with assumed scale height H
    if np.all( [varname in variables.keys() for varname in ['air_pressure','air_pressure_at_sea_level']] ):
        pmsl = variables['air_pressure_at_sea_level'][0,:,:] # [lev,lat,lon]
        zth = np.zeros(np.shape(p))
        for tstep in range(2):
            zth[tstep] = -(287*300/9.8)*np.log(p[tstep]/pmsl[np.newaxis,:,:])
    
        variables['zth']= zth

    # Potential temperature based on https://en.wikipedia.org/wiki/Potential_temperature
    # with gas constant R = 287.05 and specific heat capacity c_p = 1004.64
    if np.all( [varname in variables.keys() for varname in ['air_pressure','air_temperature']] ):
        theta = np.zeros(np.shape(p))
        Ta  = variables['air_temperature'][:,:] # [lat,lon] at surface
        repTa = np.repeat(Ta[np.newaxis,:,:], nz, axis=0) # repeat Ta along z axis
        for tstep in range(2):
            theta[tstep] = repTa*(1e5/p[tstep])**(287.05/1004.64)
        variables['theta'] = theta

    ## Destagger winds
    #u = np.tile(np.nan,(nz,ny,nx))
    if  np.all( [varname in variables.keys() for varname in ['x_wind','y_wind']] ):
        u1  = variables['x_wind'][:,:,:,:] #[time,levs,lat,lon1] wind speeds are on their directional grid edges
        v1  = variables['y_wind'][:,:,:,:] #[time,levs,lat1,lons]
        u = np.tile(np.nan,u1.shape) # tile repeats the nan accross nz,ny,nx dimensions
        u[:,:,:,1:] = 0.5*(u1[:,:,:,1:] + u1[:,:,:,:-1]) # interpolation of edges
        v = 0.5*(v1[:,:,1:,:] + v1[:,:,:-1,:]) # interpolation of edges
        s = np.hypot(u,v) # Speed is hypotenuse of u and v
        ## HYPOT on this 
        # S[0,:,:,0] ARE ALL -5000
        # S[1,:,:,0] ARE ALL NaN
        s[:,:,:,0] = np.NaN # set that edge to NaN
    
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
    if isinstance(fpaths, str):
        fpaths = [fpaths]
    # dictionary of data to keep
    data = {}
    
    ## first grab topography from pa vile
    topog,latt,lont = read_topog(_topog_sirivan_)
    data['topog']=topog
    
    ## read first file
    data0 = read_pcfile(fpaths[0],keepvars=keepvars)
    
    # Keep dimensions:
    dims = ['latitude','longitude',
            'time_1', # how many time steps in file 
            'time_0', # time step of file
            'model_level_number', # 140 levels
            'pseudo_level', # 5 levels for soil?
            'longitude_0','latitude_0', # staggered lat,lon
            'time_bnds', # file time bounds
            ]
    # copy dims
    for dim in dims:
        if dim in data0.keys():
            data[dim] = data0[dim]
    
    #assert np.all(data['latitude'] == latt), "topog latitude doesn't match pc file latitude"
    #assert np.all(data['longitude'] == lont), "topog longitude doesn't match pc file longitude"
    if (not np.all(data['latitude'] == latt)) or (not np.all(data['longitude'] == lont)):
        data['latt'] = latt
        data['lont'] = lont
    
    ## data without time dim
    flatdata = ['air_temperature',      # [lat, lon] K at surface
                'specific_humidity',    # [lat,lon] kg/kg
                'surface_air_pressure', # [lat,lon] Pa
                'surface_temperature',  # [lat,lon] K
                #'level_height',         # [z] metres from ground level]
                ]
    if keepvars is not None:
        flatdata = list(set(flatdata) & set(keepvars))
    
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
    if 'time_0' in data0.keys():
        hour.append(data0['time_0']) 
    else:
        hour.append(np.nanmean(data0['time_1']))
    time.append(data0['time_1'][0])
    time.append(data0['time_1'][1])
    
    ## Copy all the data we want
    for varname in flatdata:
        # add time dim to the flat data
        data[varname] = data0[varname][np.newaxis]
    ## And read the stuff with time dim normally
    for varname in fulldata:
        data[varname] = data0[varname]
        
    ## also copy height
    if 'level_height' in data0.keys():
        data['level_height'] = data0['level_height'][:toplev]
    
    del data0
    
    # read each file and append along time dimension
    for i,fpath in enumerate(fpaths):
        if i==0:
            continue
        datai = read_pcfile(fpath,keepvars=keepvars)
        
        for varname in flatdata:
            # add time dim to the flat data
            datai[varname] = datai[varname][np.newaxis]
        
        # Append each file along the time dimension, which has been added to the vars missing a time dim
        for varname in flatdata+fulldata:
            print("DEBUG:",varname,data[varname].shape, datai[varname].shape)
            data[varname] = np.append(data[varname], datai[varname], axis=0)
            print("DEBUG:",varname,data[varname].shape)
        
        # Append time step
        hour.append(datai['time_0']) 
        time.append(datai['time_1'][0])
        time.append(datai['time_1'][1])
        
        del datai


    
    data['hour'] = np.array(hour)
    data['time'] = np.array(time)
    # also for convenience add cloud metric
    if np.all( [varname in data.keys() for varname in ['mass_fraction_of_cloud_ice_in_air', 'mass_fraction_of_cloud_liquid_water_in_air']] ):
        data['qc'] = data['mass_fraction_of_cloud_ice_in_air'] + data['mass_fraction_of_cloud_liquid_water_in_air']


    # Finally chop the top
    for varname in fulldata:
		# Just read up to toplev 
        data[varname] = data[varname][:,:toplev,:,:]
    

    return data





def read_waroona(dtime,slv=True,ro=True,th1=True,th2=True):
    '''
        Read the converted waroona model output files
        returns list of 4 file output dictionaries
    '''
    dstamp=dtime.strftime('%Y%m%d%H')
    
    j=[slv,ro,th1,th2]
    # 4 files per datetime to be read
    keepdata = [{},{},{},{}] # keep 4 dictionaries for output
    for i,(filetype,varnames) in enumerate(_file_types_waroona_.items()):
        
        if j[i]:
            path = _files_waroona_%(dstamp,filetype)
            data=read_nc(path)
            #print("DEBUG: ",path)
            #print("DEBUG: ",data.keys())
            
            for varname in varnames:
                keepdata[i][varname] = data[varname]
            del data # delete from ram the stuff we don't want to keep
    
    # Potential temperature based on https://en.wikipedia.org/wiki/Potential_temperature
    # with gas constant R = 287.05 and specific heat capacity c_p = 1004.64
    if th1:
        k=2
        p=keepdata[k]['air_pressure']
        nt,nz,ny,nx = p.shape
        theta = np.zeros(np.shape(p))
        Ta  = keepdata[k]['air_temperature'][:,0:1,:,:] # [t,1,lat,lon] at surface
        repTa = np.repeat(Ta[:,:,:,:], 140, axis=1) # repeat Ta along z axis
        for tstep in range(nt):
            theta[tstep] = repTa*(1e5/p[tstep])**(287.05/1004.64)
        keepdata[k]['theta'] = theta

    ## Destagger winds
    #u = np.tile(np.nan,(nz,ny,nx))
    if  ro:
        k=1
        u1  = keepdata[k]['x_wind'][:,:,:,:] #[time,levs,lat,lon1] wind speeds are on their directional grid edges
        v1  = keepdata[k]['y_wind'][:,:,:,:] #[time,levs,lat1,lons]
        u = np.tile(np.nan,u1.shape) # tile repeats the nan accross nz,ny,nx dimensions
        u[:,:,:,1:] = 0.5*(u1[:,:,:,1:] + u1[:,:,:,:-1]) # interpolation of edges
        v = 0.5*(v1[:,:,1:,:] + v1[:,:,:-1,:]) # interpolation of edges
        s = np.hypot(u,v) # Speed is hypotenuse of u and v
        ## HYPOT on this 
        # S[0,:,:,0] ARE ALL -5000
        # S[1,:,:,0] ARE ALL NaN
        s[:,:,:,0] = np.NaN # set that edge to NaN
        # could also jsut set them to the adjacent edge
        #s[:,:,:,0] = s[:,:,:,1]
        assert np.all(not np.isnan(s[:,:,:,1:])), "Some nans are left in the wind_speed calculation"
        assert np.all(not s==-5000), "Some -5000 values remain in the wind_speed calculation"
        
        keepdata[k]['x_wind_destaggered'] = u
        keepdata[k]['y_wind_destaggered'] = v
        keepdata[k]['wind_speed'] = s
        
    return keepdata

def read_waroona_old(fpaths,toplev=80,keepvars=None,old=True):
    '''
    Read converted waroona output
    '''
    # converted output used to match sirivan output exactly:
    data= read_sirivan(fpaths,toplev,keepvars)
    data['topog'], latt, lont = read_topog(_topog_waroona_)
    del data['latt']
    del data['lont'] 
     
    return data
    
    
def read_topog(pafile='data/waroona/umnsaa_pa2016010515.nc'):
    '''
    Read topography and lat,lon from pa file
    '''
    with Dataset(pafile,'r') as ncfile:
        topog = ncfile.variables['surface_altitude'][:,:]
        latt = ncfile.variables['latitude' ][:]
        lont = ncfile.variables['longitude'][:]
    return topog,latt,lont
