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
from datetime import datetime
from utilities import utils, fio, constants, plotting


### GLOBALS
_sn_='localtests'



### CODE
from metrics import metric_file_path


mr="waroona_run3"
extent="waroonaz"
fpath=metric_file_path(mr,extent)

ds = xr.open_dataset(fpath)
print("INFO: reading/plotting ",fpath)
print("DEBUG:",ds)
plt.subplot(3,1,1)
ds.s_mean.isel(level=0).plot()
plt.subplot(3,1,2)
ds.firespeed_mean.plot.line("b--^")
#plt.subplot(3,1,5)
ds.firespeed_nonzero_mean.plot(color='r',linewidth=5)
print(ds.firespeed_nonzero_mean)
plt.subplot(3,1,3)
ds.sensibleheat_mean.plot.line("b--^")
print(ds.sensibleheat_nonzero_mean)
#plt.subplot(3,2,6)
ds.sensibleheat_nonzero_mean.plot(color='r',linewidth=5)
fio.save_fig_to_path("test.png",plt)