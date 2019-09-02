#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""

import matplotlib as mpl
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
from netCDF4 import Dataset

from datetime import datetime, timedelta
import timeit

import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt
import matplotlib.patches as mpatches

from utilities import utils,fio,plotting

# Try reading file using iris
import iris
import iris.quickplot as qplt


#### SETUP PLOTTING:
plotting.init_plots()

### DATETIMES FOR UM AND FIRE OUTPUT
um_dtimes = np.array([datetime(2016,1,5,15,10) + timedelta(hours=x) for x in range(24)])

dtime0 = um_dtimes[0]
dtimeE = um_dtimes[-1]
um_hour=datetime(dtime0.year,dtime0.month,dtime0.day,dtime0.hour)
ff_dtimes = np.array([um_hour + timedelta(hours=x/60.) for x in range(10,61,10)])


##### TEST CODE HERE

# get plot extent, and transect
extentname='waroona'
extent = plotting._extents_[extentname]
transect=3
start,end = plotting._transects_["%s%d"%(extentname,transect)]

## Read topography, and winds
tstep = 0
slv,ro1,th1,th2 = fio.read_waroona(dtime0, extent=extent, add_winds=True)
z, = ro1.extract('z_ro')
topog, = slv.extract('topog')
u,v,s = ro1.extract(['u','v','s'])
qc, = th2.extract('qc')
w, = th1.extract('upward_air_velocity')


# Get sign and magnitude of u,v along x axis (transect vector)
# u is along lons, v is along lats, start,end are [lat,lon]
ui = u[tstep].data
vi = v[tstep].data
lats,lons = u.coord('latitude').points,u.coord('longitude').points
ut = utils.cross_section(ui,lats,lons,start,end,npoints=80)
vt = utils.cross_section(vi,lats,lons,start,end,npoints=80)
Xprojt = (vi*(end[1]-start[1]) + ui*(end[0]-start[0])) / np.sqrt((end[0]-start[0])**2+(end[1]-start[1])**2)
Xproj0 = Xprojt * (end[0]-start[0])
Xproj1 = Xprojt * (end[1]-start[1])
X=np.sqrt(Xproj0**2 + Xproj1**2)  # horizontal magnitude along slice dim
negs=Xproj0+Xproj1 < 0.0 # where is it negative
X[negs] = -X[negs]
tranX = utils.cross_section(X,lats,lons,start=start,end=end,npoints=80)


plt.figure(figsize=[7,10])
plt.subplot(3,1,1)

# top panel is topography
plotting.map_topography(extent,topog,lat,lon)
plt.title('Topography, winds')

# Add fire front contour
with warnings.catch_warnings():
    # ignore warning when there are no fires:
    warnings.simplefilter('ignore')
    plt.contour(lon,lat,np.transpose(ff),np.array([0]), 
                colors='red')

# start to end x=[lon0,lon1], y=[lat0, lat1]
plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
         linewidth=2, 
         marker='X', markersize=7,markerfacecolor='white')

# add nearby towns
textcolor='k'
if extentname == 'waroona':
    plotting.map_add_locations(['waroona','yarloop'], 
                               text=['Waroona', 'Yarloop'], 
                               textcolor=textcolor)
    # add fire ignition
    plotting.map_add_locations(['fire_waroona'],
                               text = ['Fire ignition'], 
                               color='r', marker='*', 
                               textcolor=textcolor)
    # add pyroCB
else:
    plotting.map_add_locations(['sirivan','uarbry'], 
                               text=['Sir Ivan','Uarbry'],
                               dx=[-.02,.05], dy =[-.015,-.03],
                               textcolor=textcolor)
    # add fire ignition
    plotting.map_add_locations(['fire_sirivan'],
                               text = ['Fire ignition'], dx=.05,
                               color='r', marker='*', 
                               textcolor=textcolor)
    # add pyroCB


# Add vectors for winds
# just surface, and one every 10 points to reduce density
skip = (slice(None,None,vectorskip),slice(None,None,vectorskip))
## colour the arrows
# map wind speed to colour map domain [0, 1]
norm = col.Normalize()
norm.autoscale(s[skip])
plt.get_cmap(plotting._cmaps_['windspeed'])
plt.quiver(lon[skip[1]],lat[skip[0]],u[0][skip],v[0][skip], scale=quiverscale)
           #color=cmap(norm(s[skip])), 


## Second row is transect plots
plt.subplot(3,1,2)
wslice, xslice, zslice  = plotting.transect_w(w,z, lat, lon,start,end,
                                              npoints=100,topog=topog,
                                              lines=None)
plt.ylabel('height (m)')
## Add contour where clouds occur
qcslice = utils.cross_section(qc,lat,lon,start,end, npoints=100)
with warnings.catch_warnings():
    # ignore warning when there are no clouds:
    warnings.simplefilter('ignore')
plt.contour(xslice,zslice,qcslice,np.array([cloud_threshold]),colors='teal')

ax3 = plt.subplot(3,1,3)
trs, trx, trz = plotting.transect_s(s,z,lat,lon,start,end,topog=topog)
xticks,xlabels = plotting.transect_ticks_labels(start,end)
plt.xticks(xticks[0::2],xlabels[0::2]) # show transect start and end

# Annotate max wind speed
# only care about lower troposphere
# 60 levels is about 2800m, 80 levels is about 4700m, 70 levels : 3700m
upto = 70 
mloc = np.unravel_index(np.argmax(trs[:upto,:],axis=None),trs[:upto,:].shape)
note="max = %5.1f m/s\n  (at %4.0f m) "%(trs[:upto,:][mloc], trz[:upto,:][mloc])
trans = ax3.get_xaxis_transform() # x in data untis, y in axes fraction
ax3.annotate(note, fontsize=15,
             xy=(0.33, -0.2 ), xycoords=trans)

# Save figure into animation folder with numeric identifier
plt.suptitle(stitle)
print("INFO: Saving figure:",pname)
plt.savefig(pname)



