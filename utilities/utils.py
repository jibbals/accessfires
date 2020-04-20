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

def distance_between_points(latlon0,latlon1):
    """
    return distance between lat0,lon0 and lat1,lon1
    
    calculated using haversine formula, shortest path on a great-circle
     - see https://www.movable-type.co.uk/scripts/latlong.html
    """
    R = 6371e3 # metres (earth radius)
    lat0, lon0 = latlon0
    lat1, lon1 = latlon1
    latr0 = np.deg2rad(lat0)
    latr1 = np.deg2rad(lat1)
    dlatr = np.deg2rad(lat1-lat0)
    dlonr = np.deg2rad(lon1-lon0)
    a = np.sin(dlatr/2.0)**2 + np.cos(latr0)*np.cos(latr1)*(np.sin(dlonr/2.0)**2)
    c = 2*np.arctan2(np.sqrt(a),np.sqrt(1-a))
    return R*c

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

def cube_to_xyz(cube,
                ztop=-1):
    """
    take iris cube [lev, lat, lon]
    pull out the data and reshape it to [lon,lat,lev]
    """
    assert len(cube.shape)==3, "cube is not 3-D"
    data = cube[:ztop,:,:].data.data
    # data is now a level, lat, lon array... change to lon, lat, lev
    xyz = np.moveaxis(np.moveaxis(data,0,2),0,1)
    return xyz

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

def dates_from_iris(timedim, remove_seconds=True):
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
        for i,d in enumerate(dt):
            iday = d.day
            ihour = d.hour
            iminute = d.minute+int(d.second>30)
            if iminute == 60:
                ihour+=1
                iminute=0
            if ihour == 24:
                iday+=1
                ihour=0
            dtm[i] = datetime(d.year,d.month, iday, ihour, iminute)
    return dtm

def firepower_from_cube(shcube):
    """
    calculate and return firepower over time in GWatts
    
    Inputs
    ======
        shcube: sensible heat flux in Watts/m2
    """
    
    lon,lat = shcube.coord('longitude'), shcube.coord('latitude')
    ### get areas in m2
    # Add boundaries to grid
    if lat.bounds is None:
        lat.guess_bounds()
    if lon.bounds is None:
        lon.guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(shcube)

    firepower = shcube.data.data * grid_areas # W/m2 * m2
    return firepower/1e9 # Watts to Gigawatts

def lat_lon_index(lat,lon,lats,lons):
    ''' lat,lon index from lats,lons    '''
    with np.errstate(invalid='ignore'):
        latind=(np.abs(lats-lat)).argmin()
        lonind=(np.abs(lons-lon)).argmin()
    return latind,lonind

def nearest_date_index(date, dates, allowed_seconds=120):
    """
    Return date index that is within allowewd_seconds of date
    """
    # Need to use total seconds, close negative time deltas have -1 day + ~80k seconds
    secs_diff = np.abs([ tdelta.total_seconds() for tdelta in (np.array(dates) - date)])
    ind = np.argmin(secs_diff)
    assert secs_diff[ind] <= allowed_seconds, "%s not within %d seconds of %s ... %s. \n "%(date.strftime("%Y%m%d-%H:%M"), allowed_seconds, dates[0].strftime("%Y%m%d-%H:%M"), dates[-1].strftime("%Y%m%d-%H:%M")) + str(dates[:])
    return ind

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
    # maybe check for units
    if np.max(p) < 50000:
        print("WARNING: Potential temperature assumes pressure is in Pascals, highest input pressure is only ", np.max(p))
    if np.min(T) < 50:
        print("WARNING: Potential temperature assumes Temperature is in Kelvin, lowest input temp is only ", np.min(T))
    return T*(1e5/p)**(287.05/1004.64)
    
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

def vorticity(u,v,lats,lons):
    """
    
    ARGUMENTS:
        u = longitudinal wind [lats,lons] (m/s)
        v = latitudinal wind [lats,lons] (m/s)
        lats (deg)
        lons (deg)
        
    RETURNS: zeta is transformed from m/s/deg to 1/s
        zeta, OW_norm, OWZ
        
    NOTES:
        derivatives calculated using numpy gradient function
            >>> f = np.array([1, 2, 4, 7, 11, 16], dtype=float)
            >>> np.gradient(f)
            array([1. , 1.5, 2.5, 3.5, 4.5, 5. ])
            ## array locations can be inserted , essentially used as devision terms
            ## gradient at each point is mean of gradient on either side
            ## eg: 2 dim array
            np.gradient(np.array([[1, 2, 6], [3, 4, 5]], dtype=float))
                [array([[ 2.,  2., -1.], [ 2.,  2., -1.]]), # first output is along rows
                array([[1. , 2.5, 4. ], [1. , 1. , 1. ]])] # second output is along columns
            # This method has a dummy test in tests.py -> vorticity_test()
        Hi Jesse,
        Vorticity = zeta = v_x - u_y (where v_x = dv/dx. u_y = du/dy).
        Shearing deformation = F = v_x + U_y
        Stretching deformation = E = u_x - v_y
        OW = zeta^2 - (E^2 + F^2)        Here "^2" means squared
        OW_norm = OW/zeta^2
        OWZ = OW/zeta
        Try plotting zeta, OW_norm and OWZ.  They should all be insightful.
        Cheers,
        Kevin.
        
    """
    lat_deg_per_metre = 1/111.32e3 # 111.32km per degree
    lat_mean = np.mean(lats)
    lon_deg_per_metre = lat_deg_per_metre * np.cos(np.deg2rad(lat_mean))
    
    mlats = lats / lat_deg_per_metre # convert lats into metres
    mlons = lons / lon_deg_per_metre # convert lons into metres
    
    # u[lat,lon]
    u_lat, u_lon = np.gradient(u,mlats,mlons)
    v_lat, v_lon = np.gradient(v,mlats,mlons)
    # u is left to right (longitudinal wind)
    # v is south to north (latitudinal wind)
    zeta = v_lon - u_lat
    F = v_lon + u_lat
    E = u_lon - v_lat
    OW = zeta**2 - (E**2 + F**2)
    OW_norm = OW/(zeta**2)
    OWZ = OW/zeta
    return zeta, OW_norm, OWZ