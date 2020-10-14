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
    # FFDI is single leve, based on 10m winds and surface T,RH
    dimarraypairs["FFDI_mean"]=("time",np.zeros([ntimes])+np.NaN)
    dimarraypairs["FFDI_5ns"]=(("time","pctl"),np.zeros([ntimes,5])+np.NaN)

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
                   add_topog=True, add_z=False, add_winds=True, add_theta=False,
                   add_dewpoint=False, add_RH=True, HSkip=None, constraints=None)
    
    model_heights = utils.height_from_iris(cubes.extract("air_temperature")[0])
    
    return_dict = {}
    for varname in [
        "air_temperature", # temperature
        "air_pressure", # pressure
        "wind_direction", # wind direction
        "s", # wind speed
        "relative_humidity", # rel humidity
        ]:
        print("INFO: collating ",varname)

        cube0 = cubes.extract(varname)[0]
        # interp to level heights
        coord_names=[coord.name() for coord in cube0.coords()]
        height_coord_name = None
        if 'level_height' in coord_names:
            height_coord_name='level_height'
        elif 'atmosphere_hybrid_height_coordinate' in coord_names:
            height_coord_name='atmosphere_hybrid_height_coordinate'
        
        if height_coord_name is None:
            # get closest height indices
            hinds=np.zeros(len(heights)).astype(int)
            for i,wanted_height in enumerate(heights):
                hind = np.argmin(np.abs(model_heights-wanted_height))
                hinds[i]=hind
            #print("DEBUG: closest indices:",hinds)
            # subset to height indices
            cube  = cube0[:,hinds,:,:]
        else:
            cube = iris.util.squeeze(
                cube0.interpolate(
                    [(height_coord_name, heights)],
                    iris.analysis.Linear()
                    )
                )
        if varname=="s":
            WS_10m=np.squeeze(cube[:,1,:,:].data) # 10m wind speed
        
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
    # Surface RH as %
    RH=100*np.squeeze(cubes.extract('relative_humidity')[0][:,0,:,:].data)
    # Surface T in Celcius
    T=np.squeeze(cubes.extract('air_temperature')[0][:,0,:,:].data) - 273.15
    # WS in km/h
    FFDI=utils.FFDI(DF,RH,T,WS_10m*3.6)
    #print("DEBUG: FFDI params:")
    #print(np.mean(RH))
    #print(np.mean(T))
    #print(np.mean(WS_10m*3.6))
    #print(np.mean(FFDI))
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
            print("DEBUG: updating ",varname)
            # time, level, [pctl]
            ds[varname][insertinds] = arrays[varname]
        
        print("INFO: saving to temporary file path ",fpath_tmp)
        ds.to_netcdf(path=fpath_tmp,mode="w")# write to temporary path
    print("INFO: overwriting file path ",fpath)
    shutil.copy(fpath_tmp,fpath) # overwrite original with updated
    


mr = "sirivan_run5"
extent="sirivanz"
fpath=metric_file_path(mr, extent)
fpath = make_empty_metrics_file(mr,extentname=extent)
for hour in range(24):
    add_to_metrics_file(mr, hour=1, extentname=extent)

with xr.open_dataset(fpath) as ds:
    
    print(" ===                    === ")
    print(" === Reading file       === ")
    print(" ===                    === ")
    print(ds.keys())
    for key in ['firespeed','s','firespeed_nonzero']:
        plt.plot(ds[key], label=key,)
    plt.save_fig("test_metric.png")

        
