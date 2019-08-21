#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 11:39:09 2019

@author: jesse
"""

import numpy as np
from datetime import datetime,timedelta

# interpolation package
from scipy import interpolate


# more print statements for testing
__VERBOSE__=False


def cross_section(data,lats,lons,start,end,npoints=None):
  '''
    interpolate along cross section
    inputs: data = array[[z],lats,lons]
            start = [lat0,lon0]
            end = [lat1,lon1]
            nx = how many points along horizontal?
              If this is None, use resolution of grid
  '''
  
  lat1,lon1=start
  lat2,lon2=end
  
  # add z axis if there is none
  if len(data.shape) < 3:
    data=data[np.newaxis,:,:]
  nz = data.shape[0]

  # base interp points on grid size
  if npoints is None:
    dy = lats[1] - lats[0]
    dx = lons[1] - lons[0]
    #print (dy, dx)
    dgrid = np.hypot(dx,dy)
    #print(dgrid)
    dline = np.sqrt((lon2-lon1)**2 + (lat2-lat1)**2) # dist from start to end
    npoints = int(np.ceil(dline/dgrid))
  
  # Define grid for horizontal interpolation. x increases from 0 to 1 along the
  # desired line segment
  slicex = np.linspace(0.0,1.0,npoints)
  slicelon = lon1 + (lon2-lon1)*slicex
  slicelat = lat1 + (lat2-lat1)*slicex
  

  # Interpolate data along slice (in model height coordinate). Note that the
  # transpose is needed because RectBivariateSpline assumes the axis order (x,y)
  slicedata = np.tile(np.nan, [nz, npoints])
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
  
  return np.squeeze(slicedata)


def date_from_gregorian(greg, d0=datetime(1970,1,1,0,0,0)):
    '''
        gregorian = "hours since 19700101 00:00:00"
        Returns nparray of datetimes
    '''
    greg=np.array(greg)
    if greg.ndim==0:
        return np.array( [d0+timedelta(seconds=int(greg*3600)),])
    return np.array([d0+timedelta(seconds=int(hr*3600)) for hr in greg])


def lat_lon_index(lat,lon,lats,lons):
    ''' lat,lon index from lats,lons    '''
    with np.errstate(invalid='ignore'):
        latind=(np.abs(lats-lat)).argmin()
        lonind=(np.abs(lons-lon)).argmin()
    return latind,lonind


def relative_humidity_from_specific(qair, temp, press = 1013.25):
    '''
    modified from https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
    qair specific humidity, dimensionless (e.g. kg/kg) ratio of water mass / total air mass
    temp degrees K
    press pressure in mb
    
    return rh relative humidity, ratio of actual water mixing ratio to saturation mixing ratio
    Author David LeBauer
    '''
    tempC= temp-273.16
    es =  6.112 * np.exp((17.67 * tempC)/(tempC + 243.5))
    e  = qair * press / (0.378 * qair + 0.622)
    rh = e / es
    rh[rh > 1] = 1
    rh[rh < 0] = 0
    return(rh)
