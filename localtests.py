#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""

#import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker, colors, patches
import numpy as np
import iris
from scipy.stats import gaussian_kde, cumfreq

from datetime import datetime, timedelta

from utilities import fio, plotting, utils

_sn_ = 'localtests'

import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from owslib.wmts import WebMapTileService

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy


####
# Check hr run output
####
run_name = 'sirivan_run2_hr'
run_name2 = 'sirivan_run2'
run_hour = datetime(2017,2,11,21)

ff, = fio.read_fire(run_name)
print(ff)


sir2hrsub = fio.read_model_run(run_name,run_hour,HSkip=3, add_winds=True)
print(sir2hrsub[0].summary(shorten=True))


