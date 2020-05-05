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
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import iris # so we can add level number constraint

# read tiff file
from osgeo import gdal, osr
import cartopy.crs as ccrs

from utilities import utils, plotting, fio

model_run='waroona_run1'
extentname=model_run.split('_')[0]
HSkip=None
extent=plotting._extents_[extentname]
localtime = False

## read all the fire data, and lats/lons
sh, = fio.read_fire(model_run=model_run, dtimes=None,
                    extent=extent, firefront=False, sensibleheat=True)
lons,lats = sh.coord('longitude').points, sh.coord('latitude').points
assert sh is not None, "Missing sensible heat file"
# datetimes list, convert to local time if desired
ftimes = utils.dates_from_iris(sh)
if localtime:
    offset = timedelta(hours=8)
    if extent[0] > 130:
        offset = timedelta(hours=10)
    ftimes = np.array([ft + offset for ft in ftimes ])

# spatial sum over time of sensible heat gives firepower
firepower = np.sum(utils.firepower_from_cube(sh), axis=(1,2))

# remove zeros:
prefire = np.isclose(np.cumsum(firepower), 0)
firepower[prefire] = np.NaN

        
## Read PFT
# calculated using kevin's code
lat,lon = plotting._latlons_["fire_%s_upwind"%extentname]
pft, ptimes, _, _ = fio.read_pft(model_run,lats=lat,lons=lon)
#pft = np.nanmean(pft,axis=(1,2)) # average spatially
if localtime:
    ptimes = np.array([pt + offset for pt in ptimes ])

## Read Temperature:
# just want surface for most of these, lets just take bottom two rows:
# surf_constraint = iris.Constraint(model_level_number = lambda cell: cell<=2)
#cubes=fio.read_model_run(model_run, extent=extent, constraints=surf_constraint, 
#                         add_RH=True)

# READ EVERYTHING, SUBSET, CALC DESIRED METRICS, PLOT METRICS
cubes = fio.read_model_run(model_run, extent=extent)
ctimes = utils.dates_from_iris(cubes[0])

## get temperature, RH
q,T = cubes.extract(['specific_humidity','air_temperature'])
# compute RH from specific and T in kelvin
T.convert_units('K')
# just want surface for Temp and RH
q = q[:,0,:,:]
T = T[:,0,:,:].data.data
RH = utils.relative_humidity_from_specific(q.data.data, T)
RH = np.mean(RH, axis=(1,2))
T = np.mean(T, axis=(1,2))

## get wind speed/ wind dir
u, v = cubes.extract(['x_wind','y_wind'])
u=u[:,0,:,:]
v=v[:,0,:,:]
# DESTAGGER u and v using iris interpolate
u = u.interpolate([('longitude',v.coord('longitude').points)],
                    iris.analysis.Linear())
v = v.interpolate([('latitude',u.coord('latitude').points)],
                    iris.analysis.Linear())
ws = np.mean(utils.wind_speed_from_uv_cubes(u,v).data.data, axis=(1,2))
# mean wind direction based on mean u,v
u_mean, v_mean = np.mean(u.data.data,axis=(1,2)) , np.mean(v.data.data,axis=(1,2))
#wd = utils.wind_dir_from_uv(u_mean,v_mean)

#### PLOTTING PART

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

# figure and multiple y axes for time series
fig, axes = plt.subplots(2,1,figsize=[12,8], sharex=True)
ax_fp = axes[0] # firepower axis
color_fp = "orange"
ax_T = ax_fp.twinx() # Temperature
color_T = "red"
ax_RH = ax_fp.twinx() # Rel Humid
color_RH = "magenta"
ax_ws = axes[1] # Wind speed
color_ws = 'chartreuse'

# offsets for extra y axes, and make just the one spine visible
ax_RH.spines["right"].set_position(("axes", 1.07))
make_patch_spines_invisible(ax_RH)
ax_RH.spines["right"].set_visible(True)

## Plot firepower, and PFT on one axis
plt.sca(ax_fp)
line_pft, = plt.plot_date(ptimes, pft, '--', 
                         color=color_fp, label='PFT', alpha=0.6)
line_fp, = plt.plot_date(ftimes, firepower, '-',
                        color=color_fp, label='firepower')
## plot temperature on right axis
plt.sca(ax_T)
line_T, = plt.plot_date(ctimes, T-273.15, '-',
                       color=color_T, label="Temperature",)
## plot RH on further right axis
plt.sca(ax_RH)
line_RH, = plt.plot_date(ctimes, RH, '-',
                         color=color_RH, label="RH",)

## Wind speed, wind speed stdev, wind dir quiver
plt.sca(ax_ws)
line_ws, = plt.plot_date(ctimes, ws, '*',
                         color=color_ws, label="wind speed",)
wdnorm = np.sqrt(u_mean**2 + v_mean**2)
plt.quiver(ctimes, ws, u_mean/wdnorm, v_mean/wdnorm, 
           pivot='mid', 
           alpha=.6, 
           headwidth=2, headlength=2, headaxislength=2)


# top row:
ax_fp.set_ylabel("W/m2")
ax_fp.yaxis.label.set_color(color_fp)
ax_T.set_ylabel("Temperature (C)")
ax_T.yaxis.label.set_color(color_T)
ax_RH.set_ylabel("RH (frac)")
ax_RH.yaxis.label.set_color(color_RH)

# bottom row:
ax_ws.set_ylabel("wind speed (m/s)")
ax_ws.yaxis.label.set_color(color_ws)
ax_ws.set_xlabel("Time")

## Plot periphery
axes[0].set_title("Surface summary over time")
lines = [line_pft, line_fp, line_T, line_RH]
labels = [l.get_label() for l in lines]
ax_fp.legend(lines, labels, loc='best')
#plt.suptitle(model_run,y=0.9)

lines2 = [line_ws,]
labels2 = [l.get_label() for l in lines2]
ax_ws.legend(lines2,labels2, loc='best')

# format x-ticks date times
ax_fp.xaxis.set_major_locator(mdates.HourLocator(interval=3))
ax_fp.xaxis.set_minor_locator(mdates.HourLocator())
ax_fp.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
fig.tight_layout(rect=[0,0,.99,1]) # left, bottom, right, top
fio.save_fig("test", "localtests", 'timeseries.png',plt=plt)
