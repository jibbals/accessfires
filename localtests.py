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

import iris

cubes = fio.read_model_run('waroona_old', extent=plotting._extents_['waroonaz'],
                           fdtime=[datetime(2016,1,5,15),datetime(2016,1,5,16)],
                           add_theta=True, add_winds=True, add_RH=True, 
                           add_z=True, add_topog=True)
print(cubes)
#u, = cubes.extract('u')
#time = u.coord('time')
#sap, = cubes.extract('surface_air_pressure')
#sap0 = sap.interpolate([('time',time.points)],
#                           iris.analysis.Linear())


#sap = saps.merge_cube()
#for cube in saps:
#    cubes.remove(cube)
