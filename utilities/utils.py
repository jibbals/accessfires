#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 11:39:09 2019

@author: jesse
"""

import numpy as np

# interpolation package
from scipy import interpolate


# more print statements for testing
__VERBOSE__=True


def cross_section(data,start,end,lats,lons,npoints=None):
  '''
    interpolate along cross section
    inputs: data = array[[z],lats,lons]
            start = [lat0,lon0]
            end = [lat1,lon1]
            nx = how many points along horizontal?
              If this is None, use resolution of grid
  '''
  
  lon1,lat1=start
  lon2,lat2=end
  
  # add z axis if there is none
  if len(data.shape) < 3:
    data=data[np.newaxis,:,:]
  nz = data.shape[0]

  # base interp points on grid size
  if npoints is None:
    dy = lats[1] - lats[0]
    dx = lons[1] - lons[0]
    dgrid = np.hypot(dx,dy)
    dline = np.sqrt((lon2-lon1)**2 + (lat2-lat1)**2) # dist from start to end
    npoints = np.ceil(dgrid/dline)
  
  # Define grid for horizontal interpolation. x increases from 0 to 1 along the
  # desired line segment
  slicex = np.linspace(0.0,1.0,npoints)
  slicelon = lon1 + (lon2-lon1)*slicex
  slicelat = lat1 + (lat2-lat1)*slicex
  
  ## Physical height along slice, model height coordinates. Note that the
  ## transpose is needed because RectBivariateSpline assumes the axis order (x,y)
  #slicez = np.tile(np.nan, [nz,nx])
  #for k in range(0,nz):
  #    f = interpolate.RectBivariateSpline(lonz,latz,z[k,:,:].transpose())
  #    slicez[k,:] = f.ev(slicelon,slicelat)

  # Interpolate data along slice (in model height coordinate). Note that the
  # transpose is needed because RectBivariateSpline assumes the axis order (x,y)
  slicedata = np.tile(np.nan, [nz,npoints])
  for k in range(0,nz):
      f = interpolate.RectBivariateSpline(lons,lats,data[k,:,:].transpose())
      slicedata[k,:] = f.ev(slicelon,slicelat)
  
  if __VERBOSE__:
    # Reality check, start and end data points are close to data that is located at closest grid points to the start and end points
    i1 = np.argmin(np.abs(lons-lon1))
    i2 = np.argmin(np.abs(lons-lon2))        
    j1 = np.argmin(np.abs(lats-lat1))        
    j2 = np.argmin(np.abs(lats-lat2))
    print('VERBOSE:utils.cross_section interpolation start and end points')
    print('  Nearest neighbour vs interp data:')
    print('    {:9.2f} {:9.2f}'.format(slicedata[0,0 ],data[0,j1,i1]))
    print('    {:9.2f} {:9.2f}'.format(slicedata[0,-1],data[0,j2,i2]))
  
  return slicedata

