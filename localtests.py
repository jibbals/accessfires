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

import metrics

mr = "sirivan_run5"
extname="sirivanz"
fpath = metrics.metric_file_path(mr,extname)
with xr.open_dataset(fpath) as ds:
    
    print(" ===                    === ")
    print(" === Reading file       === ")
    print(" ===                    === ")
    print(ds.keys())
    for key in ['firespeed','firespeed_nonzero']:
        plt.plot(ds[key+"_mean"], label=key,)
    plt.legend()
    #plt.savefig("test_metric.png")
    plt.show()

        
