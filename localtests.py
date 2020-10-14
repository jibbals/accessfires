#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""


### IMPORTS

import numpy as np
import xarray as xr
import pandas as pd
import os, shutil
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import iris # interpolation etc
import time # check performance

from utilities import utils, fio, constants


### GLOBALS
_sn_='localtests'



### CODE

## create netcdf 
# create fill with empty data to be filled by hour later on

def metric_file_path(mr, extentname):
    return "data/metrics/%s_%s.nc"%(mr,extentname)

def metric_file_variables(ntimes=144):
    """
        Set up variables to be used by metric files
        Where possible names match model output named variables
        "air_temperature", # temperature
        "air_pressure", # pressure
        "wind_direction", # wind direction
        "s", # wind speed
        "relative_humidity", # rel humidity
        "FFDI", # forest fire danger index
        "level_height", # level coordinate
    """
    
    #     0: air_pressure / (Pa)                 (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    #    3: air_temperature / (K)               (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    #    6: specific_humidity / (kg kg-1)       (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    #    9: upward_air_velocity / (m s-1)       (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    #      : qc / (g kg-1)                       (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    #      : z_th / (m)                         (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    #      : u / (m s-1)                    (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    #      : v / (m s-1)                    (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    #      : s / (m s-1)                    (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    #      : wind_direction / (degrees)          (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    #      : vapour_pressure / (100 Pa)     (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    #      : dewpoint_temperature / (K)     (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    #      : relative_humidity / (1)  
    
    dimarraypairs={}
    for varnames in [
            "air_temperature", # temperature
            "air_pressure", # pressure
            "wind_direction", # wind direction
            "s", # wind speed
            "relative_humidity", # rel humidity
            ]:
        # variable holding mean value
        dimarraypairs[varnames+"_mean"]=(("time","level"),np.zeros([ntimes,10])+np.NaN)
        # variable holding min, Q1, Q2, Q3, max (5 number summary)
        dimarraypairs[varnames+"_5ns"]=(("time","level","pctl"),np.zeros([ntimes,10,5])+np.NaN)
    
    # height coordinate
    dimarraypairs["level_height"]=("level",
                 np.array([0,10, 100, 500, 1000, 2000, 4000, 6000, 10000, 15000]))
    
    # some have no level dim:
    for varname in [
            "firespeed", # velocity of fire front? m/s
            "firespeed_nonzero", # with any zeros removed
            "FFDI", # FFDI is single leve, based on 10m winds and surface T,RH
            ]:
        dimarraypairs[varname+"_mean"]=("time",np.zeros([ntimes])+np.NaN)
        dimarraypairs[varname+"_5ns"]=(("time","pctl"),np.zeros([ntimes,5])+np.NaN)
    
    return dimarraypairs

def make_empty_metrics_file(mr="sirivan_run4",
                              extentname=None, 
                              **to_netcdf_args):
    """
    create metrics file, can be filled using other method
    ARGS:
        mr: model run name from fio.run_info keys
        extentname: defaults to location+"z" (eg. sirivanz for sirivan_run5)
        other arguments to be sent to xarray.DataSet.to_netcdf, defaults:
            mode="w"
            filename="data/metrics/<mr>_<extentname>.nc"
    
    Data array: variables with _mean and _5ns as follows
        "T", # temperature (C)
        "P", # pressure (hPa)
        "WD", # wind direction (meteorological degrees)
        "WS", # wind speed (m/s)
        "RH", # rel humidity (frac?)
        "FFDI", # forest fire danger index
        "firespeed", # fire speed (m/s?)
        "firespeed_nonzero", # zeros removed
        Coords:
            "time": time
            "level": model level 0, 10m, 100m, 500m, 1000m, 2000m, 4000m, 6000m, 10000m, 15000m
            "pctl": min, Q1, median, Q3, max
    
    Returns:
        path of created file
    """
    ## Set up coordinates to be used
    sirivan_10minute_dates = pd.date_range(start="2017-02-11 11:10",end="2017-02-12 11:00",periods=144)
    waroona_10minute_dates = pd.date_range(start="2016-01-05 15:10",end="2016-01-06 15:00",periods=144)
    datecoord = waroona_10minute_dates if "waroona" in mr else sirivan_10minute_dates
    coords = {
        "time": datecoord,
        "level": np.arange(10),
        "pctl":np.arange(5),
    }
    ## Set up variables to be used
    arrays=metric_file_variables(144)
    
    
    ##  Create dataset using variables and coords
    ds = xr.Dataset(
        data_vars=arrays,
        coords=coords,
        )
    
    ## Attributes to be added!!! #TODO
    
    ## Default arguments to save the netcdf file
    if extentname is None:
        extentname=mr.split("_")[0]+"z"
    #if "group" not in to_netcdf_args:
    #    to_netcdf_args["group"]=extentname
    if "mode" not in to_netcdf_args:
        to_netcdf_args["mode"]="w"
    if "path" not in to_netcdf_args:
        to_netcdf_args["path"]=metric_file_path(mr,extentname)
    
    fio.make_folder(to_netcdf_args["path"])
    print("INFO: writing ",to_netcdf_args["path"])
    ds.to_netcdf(**to_netcdf_args)
    ds.close()
    return to_netcdf_args["path"]

def make_metrics_from_model(mr,hour=0,extentname=None,):
    """
    return dict of arrays [6, 10, [5]] with _mean and _5ns as follows:
       "air_temperature", # temperature
        "air_pressure", # pressure
        "wind_direction", # wind direction
        "s", # wind speed
        "relative_humidity", # rel humidity
        "FFDI", # forest fire danger index (based on 10m winds and surface T,RH)
    Coords:
        "time": time
        "level": model level 0, 10m, 100m, 500m, 1000m, 2000m, 4000m, 6000m, 10000m, 15000m
        "pctl": min, Q1, median, Q3, max
    """
    if extentname is None:
        extentname=mr.split("_")[0]+"z"
    extent=constants.extents[extentname]
    
    dthour=fio.run_info[mr]["filedates"][hour]
    
    vars_to_make=metric_file_variables(6)
    heights=vars_to_make['level_height'][1]
    
    # Read model data for extentname
    cubes=fio.read_model_run(mr, fdtime=dthour, extent=extent,
                             add_topog=False)
    # we subset then add wind/rel humidity fields
    #print("DEBUG: read model cubes:",cubes)
    ctimes = utils.dates_from_iris(cubes.extract('air_temperature')[0])
    model_heights = utils.height_from_iris(cubes.extract("air_temperature")[0])
    
    return_dict = {}
    
    ## Get wind speed and direction (after subsetting hopefully)
    u1, v1 = cubes.extract(['x_wind','y_wind']) #staggered but shouldn't matter for timeseries
    #t0 = time.perf_counter()
    u0 = utils.interp_cube_to_altitudes(u1,heights,model_heights, closest=True)
    #t1 = time.perf_counter()
    v0 = utils.interp_cube_to_altitudes(v1,heights,model_heights,closest=True)
    #t2 = time.perf_counter()
    # destagger
    u = u0.interpolate([('longitude',v0.coord('longitude').points)],
                       iris.analysis.Linear())
    v = v0.interpolate([('latitude',u0.coord('latitude').points)],
                       iris.analysis.Linear())
    #print("DEBUG: time with interpolation: %.5f seconds"%(t1-t0))
    #print("DEBUG: time without interpolation: %.5f seconds"%(t2-t1))
    # Get wind speed cube using hypotenuse of u,v
    firespeed,u10,v10=fio.read_fire(mr, extent=extent, dtimes=ctimes,
                          firefront=False, wind=True, firespeed=True)
    s10=utils.wind_speed_from_uv_cubes(u10,v10)
    
    s0 = utils.wind_speed_from_uv_cubes(u,v)
    s=s0.data
    s[:,1,:,:] = s10.data
    cubews = iris.cube.Cube(s,
                            var_name='wind_speed',
                            dim_coords_and_dims=[(s0.coord('time'),0),
                                                 (s0.coord('model_level_number'),1),
                                                 (s0.coord('latitude'),2),
                                                 (s0.coord('longitude'),3)]
                            )
    
    # Get wind direction using arctan of y/x
    wd = utils.wind_dir_from_uv(u.data,v.data)
    wd[:,1,:,:] = utils.wind_dir_from_uv(u10.data,v10.data)
    cubewd = iris.cube.Cube(wd,
                            var_name='wind_direction',
                            units='degrees',
                            dim_coords_and_dims=[(s0.coord('time'),0),
                                                 (s0.coord('model_level_number'),1),
                                                 (s0.coord('latitude'),2),
                                                 (s0.coord('longitude'),3)])
    
    # calculate rel humidity
    q1,T1 = cubes.extract(['specific_humidity','air_temperature'])
    # compute RH from specific and T in kelvin
    q = utils.interp_cube_to_altitudes(q1,heights,model_heights, closest=True)
    T = utils.interp_cube_to_altitudes(T1,heights,model_heights, closest=True)
    RH = utils.relative_humidity_from_specific(q.data, T.data)
    iris.std_names.STD_NAMES['relative_humidity'] = {'canonical_units': '1'}
    cubeRH = iris.cube.Cube(RH, standard_name="relative_humidity",
                               var_name="RH", units="1",
                               dim_coords_and_dims=[(q.coord('time'),0),
                                                    (q.coord('model_level_number'),1),
                                                    (q.coord('latitude'),2),
                                                    (q.coord('longitude'),3)])
    
    # also wind speed at 10m for FFDI calc
    WS_10m=np.squeeze(s[:,1,:,:]) # 10m wind speed
    # Surface RH as %
    RH_surf=100*np.squeeze(RH[:,0,:,:])
    
    # also make nonzero firespeed metric
    fs_nz = np.copy(firespeed.data)
    fs_nz[fs_nz<0.001] = np.NaN
    cubefsnz = iris.cube.Cube(fs_nz, 
                               var_name="firespeed_nonzero", 
                               units="m s-1",
                               dim_coords_and_dims=[(firespeed.coord('time'),0),
                                                    (firespeed.coord('latitude'),1),
                                                    (firespeed.coord('longitude'),2)])
    
    cubes_subset={
        "wind_direction":cubewd,
        "s":cubews,
        "relative_humidity":cubeRH,
        "firespeed":firespeed,
        "firespeed_nonzero":cubefsnz,
        }
    
    for varname in ["air_temperature","air_pressure"]:
        cube0 = cubes.extract(varname)[0]
        # interp to level heights
        cube = utils.interp_cube_to_altitudes(cube0,heights,model_heights=model_heights)
        cubes_subset[varname] = cube
    
    for varname in [
        "air_temperature", # temperature
        "air_pressure", # pressure
        "wind_direction", # wind direction
        "s", # wind speed
        "relative_humidity", # rel humidity
        "firespeed", # fire speed
        "firespeed_nonzero",
        ]:
        print("INFO: collating ",varname)
        
        cube = cubes_subset[varname]
        
        ## horizontal aggregation
        cube_mean = cube.collapsed(['latitude','longitude'], 
                                   iris.analysis.MEAN)
    
        cube_5ns = cube.collapsed(['latitude','longitude'], 
                                  iris.analysis.PERCENTILE, 
                                  percent=[0,25,50,75,100])
        
        return_dict[varname+"_mean"]=cube_mean.data
        # make pctl dimension the last one
        return_dict[varname+"_5ns"]=np.moveaxis(cube_5ns.data, 0, -1)
    
    
    DF=10.0
    print("INFO: collating FFDI (assume DF = %.1f)"%DF)
    # Surface T in Celcius
    T=np.squeeze(cubes.extract('air_temperature')[0][:,0,:,:].data) - 273.15
    # WS in km/h
    FFDI=utils.FFDI(DF,RH_surf,T,WS_10m*3.6)
    
    return_dict["FFDI_mean"]=np.mean(FFDI,axis=(1,2)) # mean over latlon
    return_dict["FFDI_5ns"]=np.moveaxis(np.percentile(FFDI,[0,25,50,75,100],
                                                      axis=(1,2)),
                                        0, -1)
    return return_dict
        

def add_to_metrics_file(mr, hour=0, extentname=None,):
    """
    take one hour of model run, get 10 minutely aggregates, update metrics file
    """
    if extentname is None:
        extentname=mr.split("_")[0]+'z'
    fpath = metric_file_path(mr,extentname)
    
    ## May need to create the file
    if not os.path.isfile(fpath):
        fpath = make_empty_metrics_file(mr,extentname)
        
    # temporary path for overwriting purposes
    fpath_tmp = fpath.split(".")[0] + "_tmp.nc"
    
    ## Now we read in the model and create our timeseries
    arrays = make_metrics_from_model(mr,hour=hour,extentname=extentname)
    nsteps=np.shape(arrays["air_temperature_mean"])[0]
    if nsteps != 6:
        print("WARNING: %d timesteps in output hour!"%nsteps)
    
    # indices for inserting new data
    insertinds=slice(hour*6,hour*6+nsteps)
    
    ## Finally open the dataset, update the metrics, save to temp path, overwrite original path
    with xr.open_dataset(fpath) as ds:
        # load data so we can overwrite it
        ds.load()
        
        # move metrics into file!
        for varname in arrays.keys():
            print("INFO: updating ",varname)
            # time, level, [pctl]
            ds[varname][insertinds] = arrays[varname]
        
        print("INFO: saving to temporary file path ",fpath_tmp)
        ds.to_netcdf(path=fpath_tmp,mode="w")# write to temporary path
    print("INFO: overwriting file path ",fpath)
    shutil.copy(fpath_tmp,fpath) # overwrite original with updated
    
mr = "sirivan_run4"
extent="sirivanz"
fpath=metric_file_path(mr, extent)
fpath = make_empty_metrics_file(mr,extentname=extent)
add_to_metrics_file(mr, hour=0, extentname=extent)
with xr.open_dataset(fpath) as ds:
    
    print(" ===                    === ")
    print(" === Reading file       === ")
    print(" ===                    === ")
    print(ds.keys())
    print(ds["FFDI_mean"][:10])
    print(ds["FFDI_5ns"][:10])
    
        
