#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker, colors, patches, dates
import numpy as np


from datetime import datetime, timedelta

from utilities import fio, plotting, utils

_sn_ = 'localtests'



ff, = fio.read_fire('sirivan_run1',dtimes=[datetime(2017,2,11,21)])

lats,lons = ff.coord('latitude').points, ff.coord('longitude').points

for metric in [lats,lons]:
    print("%.5f - %.5f, mid = %.6f"%(metric[0],metric[-1], (metric[-1]+metric[0]) / 2.0))
    print("res: %.5f, count: %5d"%(metric[1]-metric[0], len(metric)))