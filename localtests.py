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
from matplotlib import patches
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import iris # so we can add level number constraint
import iris.analysis

# read tiff file
from osgeo import gdal, osr
import cartopy.crs as ccrs

from utilities import utils, plotting, fio, constants

model_run='waroona_run1'
extentname=model_run.split('_')[0]
# zoomed extent for analysis
extentnamez = extentname + 'z'
HSkip=None
extent=plotting._extents_[extentnamez]
localtime = False
proj_latlon = ccrs.PlateCarree() # lat,lon projection
qc_thresh = constants.cloud_threshold

## read all the fire data, and lats/lons
ff,sh, = fio.read_fire(model_run=model_run, dtimes=None,
                    extent=extent, firefront=True, sensibleheat=True)
ff=ff[-1] # just final firefront is needed


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

## get temperature, RH, cloud
q,T,qc = cubes.extract(['specific_humidity','air_temperature','qc'])
clats, clons = q.coord('latitude').points, q.coord('longitude').points
# compute RH from specific and T in kelvin
T.convert_units('K')
# just want surface for Temp and RH
q = q[:,0,:,:]
T = T[:,0,:,:].data.data
qc_sum = qc.collapsed('model_level_number', iris.analysis.SUM)
qc_frac = np.sum(qc_sum.data > qc_thresh, axis=(1,2))/(len(clats)*len(clons))
qc_weight = qc_sum.collapsed(['longitude','latitude'], iris.analysis.SUM)
qc_weight = qc_weight.data/np.max(qc_weight.data)
qc_heavy = qc_weight>0.75
qc_mid = (qc_weight > 0.5) * (qc_weight < 0.75)
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
ws_all = utils.wind_speed_from_uv_cubes(u,v).data.data
ws_q1 = np.quantile(ws_all, 0.25, axis=(1,2))
ws_q3 = np.quantile(ws_all, 0.75, axis=(1,2))
ws_max = np.max(ws_all, axis=(1,2))
ws = np.mean(ws_all, axis=(1,2))
# mean wind direction based on mean u,v
u_mean, v_mean = np.mean(u.data.data,axis=(1,2)) , np.mean(v.data.data,axis=(1,2))
#wd = utils.wind_dir_from_uv(u_mean,v_mean)

#######################
#### PLOTTING PART ####
#######################

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

# figure and multiple y axes for time series
fig = plt.figure(figsize=[12,12])
dx,dy = extent[1]-extent[0], extent[3]-extent[2]
extentplus = np.array(extent)+np.array([-0.3*dx,0.3*dx,-0.2*dy,0.2*dy])

# map with extent shown
fig, ax_map, proj = plotting.map_tiff_qgis(file="%s.tiff"%extentname, 
                                           fig=fig,
                                           extent=list(extentplus), 
                                           subplot_row_col_n = [3,1,1],
                                           locnames=[extentname,]
                                           )
# Add rectangle
botleft = extent[0],extent[2]
ax_map.add_patch(patches.Rectangle(xy=botleft,
                                   width=dx,
                                   height=dy,
                                   #facecolor=None,
                                   fill=False,
                                   edgecolor='blue',
                                   linewidth=2,
                                   #linestyle='-',
                                   alpha=0.6, 
                                   transform=proj_latlon))

# add fire outline
plt.sca(ax_map)
plotting.map_fire(ff.data, lats, lons)


ax_fp = plt.subplot(3,1,2) # firepower axis
color_fp = "orange"
ax_T = ax_fp.twinx() # Temperature
color_T = "red"
ax_RH = ax_fp.twinx() # Rel Humid
color_RH = "magenta"
ax_ws = plt.subplot(3,1,3, sharex=ax_fp) # Wind speed
color_ws = 'black'
color_qc = "darkblue"
color_wd = "brown"

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
plt.plot_date(ctimes, qc_frac, 'o',
              color=color_qc, 
              fillstyle='none', 
              mec=color_qc,
              mew=1,)
line_qc, = plt.plot_date(ctimes[qc_heavy], qc_frac[qc_heavy], 'o',
                         color=color_qc,
                         label="Clouds",
                         fillstyle='full')
plt.plot_date(ctimes[qc_mid], qc_frac[qc_mid], 'o',
              color=color_qc,
              fillstyle='bottom')

## Wind speed, wind speed stdev, wind dir quiver
plt.sca(ax_ws)
line_ws, = plt.plot_date(ctimes, ws, 'o',
                         color=color_ws, label="wind speed (+IQR)",)
line_wsmax, = plt.plot_date(ctimes, ws_max, '^', 
                            color=color_ws, label="max wind speed")
# add quartiles
plt.fill_between(ctimes, ws_q3, ws_q1, alpha=0.6, color='grey')
#plt.plot_date(ctimes, ws_q1, 'v', 
#              color=color_ws)
#plt.plot_date(ctimes, ws_q3, '^', 
#              color=color_ws)

# normalize windspeed for unit length quivers
wdnorm = np.sqrt(u_mean**2 + v_mean**2)
# dont show quiver at every single point
qskip = 3
plt.quiver(ctimes[::qskip], ws_q3[::qskip], 
           u_mean[::qskip]/wdnorm[::qskip], 
           v_mean[::qskip]/wdnorm[::qskip], 
           pivot='mid', 
           alpha=.9, 
           color=color_wd,
           headwidth=3, headlength=2, headaxislength=2,
           width=.004)


# top row:
ax_fp.set_ylabel("W/m2")
ax_fp.yaxis.label.set_color(color_fp)
ax_T.set_ylabel("Temperature (C)")
ax_T.yaxis.label.set_color(color_T)
ax_RH.set_ylabel("RH and rain (frac)")
ax_RH.yaxis.label.set_color(color_RH)

# bottom row:
ax_ws.set_ylabel("wind speed (m/s)")
ax_ws.yaxis.label.set_color(color_ws)
ax_ws.set_xlabel("Time")

## Plot periphery
lines = [line_pft, line_fp, line_T, line_RH, line_qc]
labels = [l.get_label() for l in lines]
ax_fp.legend(lines, labels, loc='upper left')
#plt.suptitle(model_run,y=0.9)

lines2 = [line_ws, line_wsmax]
labels2 = [l.get_label() for l in lines2]
ax_ws.legend(lines2,labels2, loc='best')

# format x-ticks date times
ax_fp.xaxis.set_major_locator(mdates.HourLocator(interval=3))
ax_fp.xaxis.set_minor_locator(mdates.HourLocator())
ax_fp.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

ax_map.set_title("Surface area-average over time")
fig.tight_layout(rect=[0,0,.99,1]) # left, bottom, right, top
fio.save_fig("test", "localtests", 'timeseries.png', plt=plt)
