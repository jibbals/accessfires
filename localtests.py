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
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from utilities import utils, fio

### GLOBALS
_sn_='localtests'



### CODE

## create netcdf 
# make sure group="/foo/bar", mode="a" to save group foo, subgroup bar, with append mode set to not destroy the file

sirivan_10minute_dates = pd.date_range(start="2017-02-11 11:10",end="2017-02-12 11:00",periods=144)
sicoords = {
    "time": sirivan_10minute_dates,
    "level": np.arange(10),
    "pctl":np.arange(5),
}

ds = xr.Dataset(
    {
         "T_mean": (("time","level"), np.zeros([144,10])+np.NaN),
         "T_5ns": (("time","level","pctl"), np.zeros([144,10,5])+np.NaN)
     },
    coords= sicoords,
    )

# save to file under sirivanz group
# mode="w" destroys any existing data
ds.to_netcdf("data/test_file.nc",group="sirivanz",mode="w")
ds.close()

# make another group in same file...
ds = xr.Dataset({
         "T_mean": (("time","level"), np.zeros([144,10])+np.NaN),
     },
    coords= sicoords,
    )
ds.to_netcdf("data/test_file.nc",group="sirivanf",mode="a")
ds.close()

# read file
print(" ===  from disk === ")
with xr.open_dataset("data/test_file.nc",group="sirivanz") as ds_sirivanz:
    print(" sirivanz group:")
    print(ds_sirivanz)
with xr.open_dataset("data/test_file.nc",group="sirivanf") as ds_sirivanf:
    print(" sirivanf group:")
    print(ds_sirivanf)


### Next: figure out how to update values for a specific time window

print(" === setting first 6 values === ")
# if file,group,variable already exist then append along the time dimension
with xr.open_dataset("data/test_file.nc",group="sirivanz") as ds_sirivanz:
    # load data so we can overwrite it
    ds_sirivanz.load()
    
    print("T_mean:",ds_sirivanz.T_mean)
    
print(" === after closing dataset: reading from file === ")
with xr.open_dataset("data/test_file.nc",group="sirivanz") as ds:
    print(ds)
# else create file,group,variable