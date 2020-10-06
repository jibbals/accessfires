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

from utilities import utils, fio, constants

### GLOBALS
_sn_='localtests'



### CODE

## create netcdf 
# create fill with empty data to be filled by hour later on

def metric_file_path(mr, extentname):
    return "data/metrics/%s_%s.nc"%(mr,extentname)

def metric_file_variables(ntimes=144):
    ### Set up variables to be used by metric files
    dimarraypairs={}
    for varnames in [
            "T", # temperature
            "P", # pressure
            "WD", # wind direction
            "WS", # wind speed
            "RH", # rel humidity
            "FFDI", # forest fire danger index
            ]:
        # variable holding mean value
        dimarraypairs[varnames+"_mean"]=(("time","level"),np.zeros([ntimes,10])+np.NaN)
        # variable holding min, Q1, Q2, Q3, max (5 number summary)
        dimarraypairs[varnames+"_5ns"]=(("time","level","pctl"),np.zeros([ntimes,10,5])+np.NaN)
    return dimarraypairs

def create_empty_metrics_file(mr="sirivan_run4",
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
    ## add coord data
    arrays["level_height"]=("level",
                            np.array([0,10, 100, 500, 1000, 2000, 4000, 6000, 10000, 15000]))
    
    
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
    
    ds.to_netcdf(**to_netcdf_args)
    ds.close()
    return to_netcdf_args["path"]

def make_metrics_from_model(mr,hour=0,extentname=None,):
    """
    return:
        dict of arrays [6, 10, [5]]:
    """
    if extentname is None:
        extentname=mr.split("_")[0]+"z"
    extent=utils._extents_[extentname]
    dthour=fio.run_info[mr]["filedates"][hour]
    insertinds=slice(hour*6,hour*6+6)
    
    # Read model data for extentname
    cubes=fio.read_model_run(mr, fdtime=dthour, extent=extent,
                   add_topog=True, add_z=False, add_winds=False, add_theta=False,
                   add_dewpoint=False, add_RH=False, HSkip=None, constraints=None)
    
    vars_to_make=metric_file_variables(6)
    for key in vars_to_make:
        print("INFO: collating %s"%key)
        

def add_to_metrics_file(mr, hour=0, extentname=None,):
    """
    take one hour of model run, get 10 minutely aggregates, update metrics file
    """
    if extentname is None:
        extentname=mr.split("_")[0]+'z'
    fpath = metric_file_path(mr,extentname)
    
    ## May need to create the file
    if not os.path.isfile(fpath):
        fpath = create_empty_metrics_file("sirivan_run4","sirivanz")
        print("INFO: created file: %s"%fpath)
    # temporary path for overwriting purposes
    fpath_tmp = fpath.split(".")[0] + "_tmp.nc"
    
    ## Now we read in the model and create our timeseries
    arrays = get_metrics(mr,extentname)
    
    ## Finally open the dataset, update the metrics, save to temp path, overwrite original path
    with xr.open_dataset(fpath) as ds:
        # load data so we can overwrite it
        ds.load()
        
        print("saved file:")
        print(ds.keys())
        print(" ===                        === ")
        print(" === setting first 2 values === ")
        print(" ===                        === ")
        ds.T_5ns[:2] = np.zeros([2,10,5])+5
        ds.T_mean[:2] = np.zeros([2,10])+5
        
        ds.to_netcdf(path=fpath_tmp,mode="w")# write to temporary path
    shutil.copy(fpath_tmp,fpath) # overwrite original with updated
    

fpath=metric_file_path("sirivan_run4", "sirivanz")
with xr.open_dataset(fpath) as ds:
    
    print(" ===                    === ")
    print(" === Reading file       === ")
    print(" ===                    === ")
    print(ds.keys())
    for key in ds:
        print(ds[key])
        # if "_mean" in key:
        #     print(ds[key].name, ds[key][:5,0])
        # elif "_5ns" in key:
        #     print(ds[key].name, ds[key][:5,0,0])