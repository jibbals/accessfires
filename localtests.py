#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""


# IMPORTS

import numpy as np
from matplotlib import pyplot as plt

from utilities import utils

#

u = np.array([[1,2,3],[2,4,6]]) # 2 rows, 3 cols
lats = [1,2] # lats = rows
lons = [1,4,5] # lons = cols
u_y, u_x = np.gradient(u,lats,lons)


