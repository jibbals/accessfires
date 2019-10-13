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

df, dfa = fio.read_AWS_wagerup()

# Check model temperature time series vs time stamp...
import iris

r1 = fio.read_model_run('waroona_run1',add_topog=False)
# get temperature
T_surf = r1.extract('air_temperature')[0][:,0,:,:]

dts = utils.dates_from_iris(T_surf) # UTC
# WAST = UTC+8
lts = [dts[i] + timedelta(hours=8) for i in range(len(dts))]
print(T_surf)

plt.close()
plt.close()
fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
for fig, ax, X, title in zip([fig1,fig2],[ax1,ax2],[dts,lts],['UTC','Local time']):
    plt.sca(ax)
    plt.plot_date(X, np.mean(T_surf.data,axis=(1,2)))
    
    # rotate and align the tick labels so they look better
    fig.autofmt_xdate()

    # use a more precise date string for the x axis locations
    #ax.fmt_xdata = dates.DateFormatter('%Y-%m-%d')
    plt.title('mean surface temperature '+title)