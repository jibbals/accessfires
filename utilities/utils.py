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

import iris

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

def nearest_date_index(date, dates, allowed_seconds=120):
    """
    Return date index that is within allowewd_seconds of date
    """
    # Need to use total seconds, close negative time deltas have -1 day + ~80k seconds
    secs_diff = np.abs([ tdelta.total_seconds() for tdelta in (np.array(dates) - date)])
    ind = np.argmin(secs_diff)
    assert secs_diff[ind] <= allowed_seconds, "%s not within %d seconds of %s ... %s. \n "%(date.strftime("%Y%m%d-%H:%M"), allowed_seconds, dates[0].strftime("%Y%m%d-%H:%M"), dates[-1].strftime("%Y%m%d-%H:%M")) + str(dates[:])
    return ind
    

def date_index(date,dates, dn=None, ignore_hours=False):
    new_date=date
    new_dn=dn
    new_dates=dates
    if ignore_hours:
        new_date=datetime(date.year,date.month,date.day)
        if dn is not None:
            new_dn=datetime(dn.year,dn.month,dn.day)
        new_dates = [datetime(d.year,d.month,d.day) for d in dates]
    whr=np.where(np.array(new_dates) == new_date) # returns (matches_array,something)
    if len(whr[0])==0:
        print ("ERROR: ",date, 'not in', new_dates[0], '...', new_dates[-1])
    elif dn is None:
        return np.array([whr[0][0]]) # We just want the match
    else:
        whrn=np.where(np.array(new_dates) == new_dn) # returns last date match
        if len(whrn[0])==0: # last date not in dataset
            print ("ERROR: ",new_dn, 'not in', new_dates[0], '...', new_dates[-1])
        return np.arange(whr[0][0],whrn[0][0]+1)

def date_from_gregorian(greg, d0=datetime(1970,1,1,0,0,0)):
    '''
        gregorian = "hours since 19700101 00:00:00"
        Returns nparray of datetimes
    '''
    greg=np.array(greg)
    if greg.ndim==0:
        return np.array( [d0+timedelta(seconds=int(greg*3600)),])
    return np.array([d0+timedelta(seconds=int(hr*3600)) for hr in greg])

def dates_from_iris(timedim, remove_seconds=False):
    '''
    input is coord('time') and grain
    or else input a cube with a 'time' dim
    output is array of datetimes
    '''
    tdim=timedim
    if isinstance(timedim, iris.cube.Cube):
        tdim  = timedim.coord('time')
        grain = str(tdim.units).split(' ')[0]
    #elif isinstance(timedim, [iris.coords.DimCoord,iris.coords.Coord]):
    elif isinstance(timedim, iris.coords.Coord):
        grain = str(timedim.units).split(' ')[0]
    
    unitformat = '%s since %%Y-%%m-%%d %%H:%%M:%%S'%grain
    d0 = datetime.strptime(str(tdim.units),unitformat)
    secmult=1
    if grain=='minutes':
        secmult=60
    elif grain=='hours':
        secmult=60*60
    elif grain=='days':
        secmult=60*60*24
    
    dt = np.array([d0 + timedelta(seconds=secs*secmult) for secs in tdim.points])
    dtm = dt
    if remove_seconds:
        dtm = np.array([datetime(d.year,d.month,d.day,d.hour,d.minute) for d in dt])
    return dtm

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

def zth_calc(pmsl,p):
    '''
    Calculate zth from pmsl and p
    pmsl : [t,lat,lon]
    p    : [t,lev,lat,lon]
    
    '''
    nt,nz,ny,nx = p.shape
    reppmsl = np.repeat(pmsl[:,np.newaxis,:,:],nz, axis=1) # repeat surface pressure along z axis
    return -(287*300/9.8)*np.log(p/reppmsl)

def potential_temperature(p,T):
    '''
    calculate theta from pressure and air temperature
    # Potential temperature based on https://en.wikipedia.org/wiki/Potential_temperature
    # with gas constant R = 287.05 and specific heat capacity c_p = 1004.64
    '''
    nt,nz,ny,nx = p.shape
    #Ta  = T[:,0:1,:,:] # [t,1,lat,lon] at surface
    #repTa = np.repeat(Ta[:,:,:,:], nz, axis=1) # repeat Ta along z axis
    #assert np.all(repTa[:,0,:,:] - repTa[:,1,:,:] == 0), "Repeated z dim is not the same"
    #return repTa*(1e5/p)**(287.05/1004.64)
    return T*(1e5/p)**(287.05/1004.64)
    #print("DEBUG: ",p.shape, theta.shape, repTa.shape, Ta.shape)
    #print("DEBUG: ",'theta' in keepdata[k].keys())
    
def destagger_winds(u1,v1):
    '''
    destagger winds from ACCESS um output
    wind speeds are on their directional grid edges
    #u1 = [time,levs,lat,lon1] 
    #v1 = [time,levs,lat1,lons]
    
    '''
    u = np.tile(np.nan,u1.shape) # tile repeats the nan accross nz,ny,nx dimensions
    u[:,:,:,1:] = 0.5*(u1[:,:,:,1:] + u1[:,:,:,:-1]) # interpolation of edges
    v = 0.5*(v1[:,:,1:,:] + v1[:,:,:-1,:]) # interpolation of edges
    
    return u,v

def wind_speed(u,v, fix=True):
    '''
    horizontal wind speed from u,v vectors
    fix sets instances of -5000 to NaN (problem from destaggering winds with one missing edge)
    '''
    
    s = np.hypot(u,v) # Speed is hypotenuse of u and v
    
    if fix and np.sum(s==-5000)>0:
        s[s==-5000] = np.NaN
        assert np.sum(np.isnan(s[:,:,:,1:]))==0, "Some nans are left in the wind_speed calculation"
        assert np.sum(s==-5000)==0, "Some -5000 values remain in the wind_speed calculation"
    
    return s
    ## HYPOT on this 
    # S[0,:,:,0] ARE ALL -5000
    # S[1,:,:,0] ARE ALL NaN
    #s[:,:,:,0] = np.NaN # set that edge to NaN
    # could also jsut set them to the adjacent edge
    #s[:,:,:,0] = s[:,:,:,1]
        
