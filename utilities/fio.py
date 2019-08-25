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
#from iris.fileformats.netcdf import load_cubes
import iris
import numpy as np
import timeit # for timing stuff
from datetime import datetime, timedelta

#from .context import utils

from glob import glob
###
## GLOBALS
###
__VERBOSE__=True

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
_file_types_waroona_ = {'slv':['specific_humidity',  # kg/kg [t,lat,lon]
                               'surface_air_pressure', # Pa 
                               'air_pressure_at_sea_level', # Pa [t,lat,lon]
                               #'surface_altitude', # [lat,lon] # topog # only occurs in first output
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
                               #'height_above_reference_ellipsoid', #[z,y,x] # new level heights? # only occurs in first outfile
                               'time', # = 6 ;
                               'level_height', # = 140 ;
                               'latitude', # = 576 ;
                               'longitude', # = 576
                               ],
                        'mdl_th1':['air_pressure', # Pa [t,z,y,x]
                                   'air_temperature', # K [t,z,y,x]
                                   #'height_above_reference_ellipsoid', # m [t,z,y,x] # only occurs in first outfile
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
    tstart = timeit.default_timer()
    ncfile =  Dataset(fpath,'r')
    print("INFO: Reading ",fpath, " ... ")
    variables = {}
    if keepvars is None:
        for vname in ncfile.variables:
            variables[vname] = ncfile.variables[vname][:]
    else:
        # keep dimensions and whatever is in keepvars
        #for vname in set(keepvars) | (set(ncfile.dimensions.keys()) & set(ncfile.variables.keys())):
        for vname in set(keepvars):
            variables[vname] = ncfile.variables[vname][:]

  
    print("INFO: finished reading %s (%6.2f minutes) "%( fpath, ( timeit.default_timer()-tstart )/60.0 ))
    return variables


def read_nc_iris(fpath, constraints=None, keepvars=None):
    '''
    Try using iris and see if it's way faster or whatever
    '''
    
    print("INFO: Reading(iris) ",fpath)
    # First read all cubes in file using constraints
    if constraints is not None:
        cubes=iris.load(fpath, constraints)
    else:
        cubes = iris.load(fpath)
    
    # If just some variables are wanted, pull them out
    if keepvars is not None:
        cubes = cubes.extract(keepvars)
    
    return cubes
    

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

def read_fire(fpath='data/waroona_fire/firefront.CSIRO_24h.20160105T1500Z.nc', 
              dtimes=None,
              cube=False):
    '''
    from dpath, read some the firefront and heat flux
    if tsteps is set to a list of datetimes, use this to make a slice of the time dimension
    '''
    # Read fire front, taking one point from each 10 minutes
    #fpath=['C:/Users/jgreensl/Desktop/data/waroona_fire/firefront.CSIRO_24h.20160105T1500Z.nc',
    #   'C:/Users/jgreensl/Desktop/data/waroona_fire/firefront.CSIRO_24h.20160106T1500Z.nc'
    #   ]
    
    print("INFO: reading fire file ",fpath, '...')
    # read the cube
    ff, = read_nc_iris(fpath)
    
    # just want one entry every 10 minutes:
    if dtimes is None:
        ff = ff[9::10]
    else:
        # We need to match the time coord in fire data to our dtimes list
        # pull datetimes from time dimension
        fft  = ff.coord('time')
        d0 = datetime.strptime(str(fft.units),'seconds since %Y-%m-%d %H:%M:00')
        # datetimes from ff
        ffdt = np.array([d0 + timedelta(seconds=secs) for secs in fft.points])
        # for each datetime in dtimes argument, find closest index in dt
        tinds = []
        for dtime in dtimes:
            tinds.append(np.argmin(abs(ffdt-dtime)))
        tinds = np.array(tinds)
        # Check that fire times are within 2 minutes of desired dtimes
        #print("DEBUG: diffs")
        #print([(ffdt[tinds][i] - dtimes[i]).seconds < 121 for i in range(len(dtimes))])
        assert np.all([(ffdt[tinds][i] - dtimes[i]).seconds < 121 for i in range(len(dtimes))]), "fire times are > 2 minutes from requested dtimes"
        
        # subset cube to desired times
        ff = ff[tinds]
        
    
    # maybe just return the iris cube
    if cube:
        return ff
    
    
    lats = ff.coord('latitude').points
    lons = ff.coord('longitude').points
    
    # here is where the data is actually read
    ffdata = ff.data
    tarr = ff.coord('time').points # seconds since d0
    
    return ffdata, tarr, lats, lons 
    
    
def read_z(fpath='data/waroona/umnsaa_2016010515_mdl_th1.nc'):
    '''
    '''
    print("INFO: reading fire file ",fpath, '...')
    # read the cube
    z, = read_nc_iris(fpath,'height_above_reference_ellipsoid')
    #height_above_reference_ellipsoid / (m) (model_level_number: 140; latitude: 576; longitude: 576)
    return z
    

def read_waroona_iris(dtime, constraints=None):
    '''
        Read the converted waroona model output files
        returns list of 4 iris cube lists:
        ========
        0: specific_humidity / (1)             (time: 6; latitude: 576; longitude: 576)
        1: surface_air_pressure / (Pa)         (time: 6; latitude: 576; longitude: 576)
        2: air_pressure_at_sea_level / (Pa)    (time: 6; latitude: 576; longitude: 576)
        3: surface_temperature / (K)           (time: 6; latitude: 576; longitude: 576)
        ========
        0: air_pressure / (Pa)                 (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        1: x_wind / (m s-1)                    (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        2: y_wind / (m s-1)                    (time: 6; model_level_number: 140; latitude: 577; longitude: 576)
        ========
        0: air_pressure / (Pa)                 (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        1: air_temperature / (K)               (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        2: specific_humidity / (kg kg-1)       (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        3: upward_air_velocity / (m s-1)       (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        ========
        0: mass_fraction_of_cloud_ice_in_air / (kg kg-1) (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        1: mass_fraction_of_cloud_liquid_water_in_air / (kg kg-1) (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
    '''
    dstamp=dtime.strftime('%Y%m%d%H')
    
    cubelists = []
    for filetype,varnames in _file_types_waroona_.items():
        path = _files_waroona_%(dstamp,filetype)
        cubelists.append(read_nc_iris(path,constraints=constraints,keepvars=varnames))

    return cubelists

def read_waroona_fire():
    '''
    '''

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
            data=read_nc(path,varnames) # only read the stuff we want!!
            #print("DEBUG: ",path)
            #print("DEBUG: ",data.keys())
            
            for varname in varnames:
                keepdata[i][varname] = data[varname]
            del data # delete from ram the stuff we don't want to keep
    

    # zth
    #
    if th1 and slv:
        pmsl = keepdata[0]['air_pressure_at_sea_level'][:,:,:] # [t,lat,lon]
        p=keepdata[2]['air_pressure'][:,:,:,:] # [t,lev,lat,lon]
        nt,nz,ny,nx = p.shape
        reppmsl = np.repeat(pmsl[:,np.newaxis,:,:],nz, axis=1) # repeat surface pressure along z axis
        zth = -(287*300/9.8)*np.log(p/reppmsl)
        keepdata[2]['zth']= zth


    # Potential temperature based on https://en.wikipedia.org/wiki/Potential_temperature
    # with gas constant R = 287.05 and specific heat capacity c_p = 1004.64
    if th1:
        k=2
        p=keepdata[k]['air_pressure']
        nt,nz,ny,nx = p.shape
        theta = np.zeros(np.shape(p))
        Ta  = keepdata[k]['air_temperature'][:,0:1,:,:] # [t,1,lat,lon] at surface
        repTa = np.repeat(Ta[:,:,:,:], nz, axis=1) # repeat Ta along z axis
        assert np.all(repTa[:,0,:,:] - repTa[:,1,:,:] == 0), "Repeated z dim is not the same"
        theta = repTa*(1e5/p)**(287.05/1004.64)
        #print("DEBUG: ",p.shape, theta.shape, repTa.shape, Ta.shape)
        #print("DEBUG: ",'theta' in keepdata[k].keys())
        if 'theta' in keepdata[k].keys():
            print("DEBUG: ", keepdata[k]['theta'].shape)
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
        assert np.sum(np.isnan(s[:,:,:,1:]))==0, "Some nans are left in the wind_speed calculation"
        assert np.sum(s==-5000)==0, "Some -5000 values remain in the wind_speed calculation"
        
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
    
    
def read_topog(pafile='data/waroona/topog.nc'):
    '''
    Read topography and lat,lon from pa file
    for new output, topography exists as 'surface_altitude' in the first slv output
    moved to topog.nc by ncks:
        ncks -v surface_altitude ..._slv.nc data/waroona/topog.nc 
    For easy reading here
    '''
    with Dataset(pafile,'r') as ncfile:
        topog = ncfile.variables['surface_altitude'][:,:]
        latt = ncfile.variables['latitude' ][:]
        lont = ncfile.variables['longitude'][:]
    return topog,latt,lont
