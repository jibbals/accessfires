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

metrics.compare_metrics(mrs=["sirivan_run4","sirivan_run5","sirivan_run6", "sirivan_run5_hr", "sirivan_run6_hr"])
