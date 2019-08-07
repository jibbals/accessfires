# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 11:22:42 2016

@author: Jeff
"""

import matplotlib.colors as col
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

def plot_xsec(lon,lat,data,lonz,latz,z,lont,latt,topog,end1,end2,clev1,clev2,ztop=None,nx=100,doplot=True,cmap='jet',norm=None,cbformat=None):
    """Plot a cross-section of some UM data on model levels
    
    lon,lat: grid that data is on
    data: 3-d data array on model levels
    z: 3-d array of physical height, same grid as data
    lont,latt: grid that topography is on
    topog: 2-d array of topography
    end1=(lon1,lat1): left-hand endpoint of section
    end2=(lon2,lat2): right-hand endpoint of section
    clev: contour levels to plot
    ztop: optional upper limit to plot
    nx: no of points used for horizontal interpolation to slice
    doplot: actually do the plot, into current axes
    
    return values: interpolated data to plot yourself
    
    Written Jeff Kepert 21 July 2016
    """

    lon1,lat1=end1
    lon2,lat2=end2
    nz = data.shape[0]

    # Define grid for horizontal interpolation. x increases from 0 to 1 along the
    # desired line segment
    slicex = np.linspace(0.0,1.0,nx)
    slicelon = lon1 + (lon2-lon1)*slicex
    slicelat = lat1 + (lat2-lat1)*slicex
    
    # Interpolate topography along slice
    f = interpolate.RectBivariateSpline(lont,latt,topog.transpose())
    slicetopog = f.ev(slicelon,slicelat)
    
    # Physical height along slice, model height coordinates. Note that the
    # transpose is needed because RectBivariateSpline assumes the axis order (x,y)
    slicez = np.tile(np.nan, [nz,nx])
    for k in range(0,nz):
        f = interpolate.RectBivariateSpline(lonz,latz,z[k,:,:].transpose())
        slicez[k,:] = f.ev(slicelon,slicelat)

    # Interpolate data along slice (in model height coordinate). Note that the
    # transpose is needed because RectBivariateSpline assumes the axis order (x,y)
    slicedata = np.tile(np.nan, [nz,nx])
    for k in range(0,nz):
        f = interpolate.RectBivariateSpline(lon,lat,data[k,:,:].transpose())
        slicedata[k,:] = f.ev(slicelon,slicelat)

    # Reality check on the interpolation, eventually remove
    i1 = np.argmin(np.abs(lont-lon1))        
    i2 = np.argmin(np.abs(lont-lon2))        
    j1 = np.argmin(np.abs(latt-lat1))        
    j2 = np.argmin(np.abs(latt-lat2))        
    print('Nearest neighbour vs interp topog:')
    print('   {:9.2f} {:9.2f}'.format(slicetopog[0 ],topog[j1,i1]))
    print('   {:9.2f} {:9.2f}'.format(slicetopog[-1],topog[j2,i2]))
    
    i1 = np.argmin(np.abs(lon-lon1))        
    i2 = np.argmin(np.abs(lon-lon2))        
    j1 = np.argmin(np.abs(lat-lat1))        
    j2 = np.argmin(np.abs(lat-lat2))        
    print('Nearest neighbour vs interp data:')
    print('   {:9.2f} {:9.2f}'.format(slicedata[0,0 ],data[0,j1,i1]))
    print('   {:9.2f} {:9.2f}'.format(slicedata[0,-1],data[0,j2,i2]))
    
    # Set up a tuple of strings for the labels. Very crude!
    xticks = (0.0,0.5,1.0)
    fmt = '{:.1f}S {:.1f}E'
    xlabs = (fmt.format(-lat1,lon1),fmt.format(-0.5*(lat1+lat2),0.5*(lon1+lon2)),fmt.format(-lat2,lon2))    
    
    if doplot:
        # Note that contourf can work with non-plaid coordinate grids provided both are 2-d
        if norm is None:
            plt.contourf(np.tile(slicex,(nz,1)),slicez,slicedata,clev1,cmap=cmap)
        else:
            plt.contourf(np.tile(slicex,(nz,1)),slicez,slicedata,clev1,cmap=cmap,norm=norm)
        if cbformat is None:
            plt.colorbar()
        else:
            plt.colorbar(format=cbformat)
        plt.contour(np.tile(slicex,(nz,1)),slicez,slicedata,clev2,colors='k')            
        plt.fill_between(slicex,slicetopog,interpolate=True,facecolor='black')
        #plt.plot(slicex,slicetopog,'k-',linewidth=2.0)
        plt.xticks(xticks,xlabs)
        if ztop != None:
            plt.ylim(0,ztop)
        plt.ylabel('Height (m)')
        plt.xlabel('Latitude and longitude')
            
    return #slicelon,slicelat,slicetopog,slicez,slicedata

def plot_xsec_spd(lonu,latu,u,lonv,latv,v,lonz,latz,z,lont,latt,topog,end1,end2,clev1,clev2,ztop=None,nx=100,doplot=True,cmap='jet'):
    """Plot a cross-section of some UM data on model levels
    
    lon,lat: grid that data is on
    data: 3-d data array on model levels
    z: 3-d array of physical height, same grid as data
    lont,latt: grid that topography is on
    topog: 2-d array of topography
    end1=(lon1,lat1): left-hand endpoint of section
    end2=(lon2,lat2): right-hand endpoint of section
    clev: contour levels to plot
    ztop: optional upper limit to plot
    nx: no of points used for horizontal interpolation to slice
    doplot: actually do the plot, into current axes
    
    return values: interpolated data to plot yourself
    
    Written Jeff Kepert 21 July 2016
    """

    lon1,lat1=end1
    lon2,lat2=end2
    nz = z.shape[0]

    # Define grid for horizontal interpolation. x increases from 0 to 1 along the
    # desired line segment
    slicex = np.linspace(0.0,1.0,nx)
    slicelon = lon1 + (lon2-lon1)*slicex
    slicelat = lat1 + (lat2-lat1)*slicex
    
    # Interpolate topography along slice
    f = interpolate.RectBivariateSpline(lont,latt,topog.transpose())
    slicetopog = f.ev(slicelon,slicelat)
    
    # Physical height along slice, model height coordinates. Note that the
    # transpose is needed because RectBivariateSpline assumes the axis order (x,y)
    slicez = np.tile(np.nan, [nz,nx])
    for k in range(0,nz):
        f = interpolate.RectBivariateSpline(lonz,latz,z[k,:,:].transpose())
        slicez[k,:] = f.ev(slicelon,slicelat)

    # Interpolate u along slice (in model height coordinate). Note that the
    # transpose is needed because RectBivariateSpline assumes the axis order (x,y)
    sliceu = np.tile(np.nan, [nz,nx])
    for k in range(0,nz):
        f = interpolate.RectBivariateSpline(lonu,latu,u[k,:,:].transpose())
        sliceu[k,:] = f.ev(slicelon,slicelat)
    slicev = np.tile(np.nan, [nz,nx])
    for k in range(0,nz):
        f = interpolate.RectBivariateSpline(lonv,latv,v[k,:,:].transpose())
        slicev[k,:] = f.ev(slicelon,slicelat)

    slicespd = np.hypot(sliceu,slicev)
    
    # Reality check on the interpolation, eventually remove
    i1 = np.argmin(np.abs(lont-lon1))        
    i2 = np.argmin(np.abs(lont-lon2))        
    j1 = np.argmin(np.abs(latt-lat1))        
    j2 = np.argmin(np.abs(latt-lat2))        
    print('Nearest neighbour vs interp topog:')
    print('   {:9.2f} {:9.2f}'.format(slicetopog[0 ],topog[j1,i1]))
    print('   {:9.2f} {:9.2f}'.format(slicetopog[-1],topog[j2,i2]))
    
    i1 = np.argmin(np.abs(lonu-lon1))        
    i2 = np.argmin(np.abs(lonu-lon2))        
    j1 = np.argmin(np.abs(latu-lat1))        
    j2 = np.argmin(np.abs(latu-lat2))        
    print('Nearest neighbour vs interp data:')
    print('   {:9.2f} {:9.2f}'.format(sliceu[0,0 ],u[0,j1,i1]))
    print('   {:9.2f} {:9.2f}'.format(sliceu[0,-1],u[0,j2,i2]))
    
    # Set up a tuple of strings for the labels. Very crude!
    xticks = (0.0,0.5,1.0)
    fmt = '{:.1f}S {:.1f}E'
    xlabs = (fmt.format(-lat1,lon1),fmt.format(-0.5*(lat1+lat2),0.5*(lon1+lon2)),fmt.format(-lat2,lon2))    
    
    if doplot:
        # Note that contourf can work with non-plaid coordinate grids provided both are 2-d
        plt.contourf(np.tile(slicex,(nz,1)),slicez,slicespd,clev1,cmap=cmap)
        plt.colorbar()
        plt.contour(np.tile(slicex,(nz,1)),slicez,slicespd,clev2,colors='k')            
        plt.fill_between(slicex,slicetopog,interpolate=True,facecolor='black')
        plt.xticks(xticks,xlabs)
        if ztop != None:
            plt.ylim(0,ztop)
        plt.ylabel('Height (m)')
        plt.xlabel('Latitude and longitude')
            
    return #slicelon,slicelat,slicetopog,slicez,slicespd

