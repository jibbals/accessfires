#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 15:14:27 2020
    Look at surfaces in 3d
        show surface level, and atmospheric stuff in a set of 3d axes
@author: jesse
"""


import matplotlib
#matplotlib.use('Agg',warn=False)

# plotting stuff
import matplotlib.pyplot as plt
import numpy as np
import warnings
from datetime import datetime

# 3d plotting
# isosurface plotting available in plotly
#import plotly.graph_objects as go
#from plotly.io import show
#import plotly.io as pio
from mayavi import mlab


# local modules
from utilities import plotting, utils, fio

###
## GLOBALS
###
_sn_ = 'threedee'

### DEMO
from numpy import pi, sin, cos, mgrid
dphi, dtheta = pi/250.0, pi/250.0
[phi,theta] = mgrid[0:pi+dphi*1.5:dphi,0:2*pi+dtheta*1.5:dtheta]
m0 = 4; m1 = 3; m2 = 2; m3 = 3; m4 = 6; m5 = 2; m6 = 6; m7 = 4;
r = sin(m0*phi)**m1 + cos(m2*phi)**m3 + sin(m4*theta)**m5 + cos(m6*theta)**m7
x = r*sin(phi)*cos(theta)
y = r*cos(phi)
z = r*sin(phi)*sin(theta)

# View it.
s = mlab.mesh(x, y, z)
mlab.show()


###
## METHODS
###


#
#
#
### Get run data
#mr='waroona_run2'
#top=10000 # km height
#hours=fio.model_outputs[mr]['filedates']
##hour=hours[0]
#hour = datetime(2016,1,6,10)
#extentname = mr.split('_')[0]
#extent = plotting._extents_[extentname]
#cubes = fio.read_model_run(mr,fdtime=hour,extent=extent, 
#                           add_theta=True, add_topog=True)
#ff, = fio.read_fire(mr,dtimes=[hour],extent=extent)
#
#print(cubes)
#
#th, qc = cubes.extract(['potential_temperature','qc'])
#topog, = cubes.extract(['surface_altitude'])
#levh  = qc.coord('level_height').points
#topind = np.sum(levh<top)
#qc0= qc[0,:topind,:,:].data.data
#th0= th[0,:topind,:,:].data.data
#levh = levh[:topind]
#
## these are level, lat, lon cubes
#lat,lon = qc.coord('latitude').points, qc.coord('longitude').points
#lev = np.arange(len(qc0[:,0,0]))
#
#
#qc0_xyz = np.moveaxis(np.moveaxis(qc0,0,2),0,1)
#th0_xyz = np.moveaxis(np.moveaxis(th0,0,2),0,1)
#topog_xy = topog.data.data.T
#ff0_xy = ff[0].data.data
#
### 3d figure:
#
#X,Y,Z = np.meshgrid(lon,lat,levh) 
#    
## X Y Z are now [lat, lon, lev] for some reason
#[X,Y,Z] = [ np.moveaxis(arr,0,1) for arr in [X,Y,Z]]
## Now they are lon, lat, lev
#
## Add topography to Z
#Zh=np.zeros(np.shape(Z))
#for i in range(np.shape(Zh)[2]):
#    Zh[:,:,i] = Z[:,:,i]+topog_xy
#
#mlab.contour3d(qc0_xyz)