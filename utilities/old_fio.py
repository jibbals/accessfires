# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 10:27:17 2019
    Old code, generally reads data without using iris and takes lots of ram/time
@author: jgreensl
"""

#import numpy as np
from netCDF4 import Dataset
#from iris.fileformats.netcdf import load_cubes
import iris
import numpy as np
import timeit # for timing stuff
from datetime import datetime, timedelta
from glob import glob


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


# Old version of waroona output:
_files_waroona_old_= sorted(glob('data/waroona/umnsaa_pc*.nc'))

def read_waroona(dtime,slv=True,ro=True,th1=True,th2=True):
    '''
        Read the converted waroona model output files
        returns list of 4 file output dictionaries
        Doesn't use iris, needs lots of RAM (10GB+)
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

def read_fire_front(fpath='data/waroona_fire/firefront.CSIRO_24h.20160105T1500Z.nc', 
                    dtimes=None):
    '''
    from dpath, read some the firefront
    if tsteps is set to a list of datetimes, use this to make a slice of the time dimension
    '''
    
    # read the cube
    ff, = read_nc_iris(fpath)
    
    if dtimes is not None:
        # We need to match the time coord in fire data to our dtimes list
        ff = subset_time_iris(ff,dtimes,seconds=True,seccheck=121)
        
    # Return the iris cube
    return ff

def read_fire_flux(fpath='data/waroona_fire/sensible_heat.CSIRO_24h.20160105T1500Z.nc', 
                    dtimes=None):
    '''
    from dpath, read some the sensible heat flux from fire output
    if tsteps is set to a list of datetimes, use this to make a slice of the time dimension
    '''
    
    print("INFO: reading fire file ",fpath, '...')
    # read the cube
    ff, = read_nc_iris(fpath)
    
    if dtimes is not None:
        ff = subset_time_iris(ff,dtimes,seconds=True,seccheck=121)
        
    # Return the iris cube
    return ff

def read_fire_old(fpath='data/waroona_fire/firefront.CSIRO_24h.20160105T1500Z.nc'):
    '''
    '''
    print("INFO: reading fire file ",fpath, '...')
    # read the cube
    ff, = read_nc_iris(fpath)
    
    # just want one entry every 10 minutes:
    ff = ff[9::10]
    
    lats = ff.coord('latitude').points
    lons = ff.coord('longitude').points
    
    # here is where the data is actually read
    ffdata = ff.data
    tarr = ff.coord('time').points # seconds since d0
    
    return ffdata, tarr, lats, lons 