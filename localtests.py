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
import iris

from datetime import datetime, timedelta

from utilities import fio, plotting, utils

_sn_ = 'localtests'

import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from owslib.wmts import WebMapTileService

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy

pft,ptimes,plats,plons = fio.read_pft()


#####
## Attempt a widget
#####
import ipywidgets as widgets
from matplotlib.widgets import Button

fig, ax = plt.subplots()
fig.subplots_adjust(bottom=0.2)

ax.hist(pft.flatten(),bins=np.arange(0,401,10))

class Index:
    bins=np.arange(0,401,10)
    data = pft
    data_min = 0
    data_max = data.shape[0]-1
    selected = 0
    def next(self, event):
        if self.selected >= self.data_max:
            self.selected = self.data_max
            ax.set_title('Last time slice reached. Cannot go forwards')
        else:
            self.selected += 1
            ax.cla()
            ax.hist(self.data[self.selected].flatten(),bins=self.bins)
            ax.set_title("time slice %d" %self.selected)

    def prev(self, event):
        if self.selected <=self.data_min:
            self.selected = 0
            ax.set_title('First time slice reached. Cannot go backwards')
        else:
            self.selected -= 1
            ax.cla()
            ax.hist(self.data[self.selected].flatten(), bins=self.bins)
            ax.set_title("time slice: %d" %self.selected)
        

callback = Index()
axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
axnext = plt.axes([0.81, 0.05, 0.1, 0.075])

bnext = Button(axnext, '>')
bnext.on_clicked(callback.next)

bprev = Button(axprev, '<')
bprev.on_clicked(callback.prev)