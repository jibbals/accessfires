#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 14:06:49 2019

  In this script is repeatable generic plotting stuff like adding a scale to a cartopy map

@author: jesse
"""

# Plotting stuff
import numpy as np
import matplotlib
from matplotlib import patheffects
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.ticker as tick

# Mapping stuff
import cartopy
import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt

# local stuff
from .context import utils

###
### GLOBALS
###

## Standard colour maps for use between figures
_cmaps_ = {}
_cmaps_['verticalvelocity'] = plt.cm.get_map('PiYG')
_cmaps_ ['windspeed'] = plt.cm.get_cmap('YlGnBu')

# Extents: EWSN
_extents_             = {}
_extents_['waroona']  = [115.775,116.2, -33.05,-32.7] # local
_extents_['waroonas'] = [112,120,-34.5,-31] # synoptic
_extents_['yarloop']  = _extents_['waroona']
_extents_['yarloops'] = _extents_['waroonas']
_latlons_             = {}
_latlons_['waroona']  = -32.8430, 115.8526 # latlon of waroona
_latlons_['yarloop']  = -32.9534, 115.9124
_latlons_['perth']    = -31.9505, 115.8605

_latlons_['fire_waroona'] = -32.89, 116.17
_latlons_['fire_sirivan'] = -30,120 # check with harvey or find in accessdev

_transects_             = {} 
_transects_['waroona1'] = [-32.75,115.75], [-32.95,116.19]
_transects_['waroona2'] = [-32.85,115.75], [-33.05,116.19] # like waroona1 but lower
_transects_['waroona3'] = [-32.65,115.75], [-32.85,116.19] # '' but higher

def init_plots():
  matplotlib.rcParams['font.size'] = 16.0
  matplotlib.rcParams["text.usetex"]      = False     # I forget what this is for, maybe allows latex labels?
  matplotlib.rcParams["legend.numpoints"] = 1         # one point for marker legends
  matplotlib.rcParams["figure.figsize"]   = (9, 7)    # Default figure size
  matplotlib.rcParams["axes.titlesize"]   = 22        # title font size
  matplotlib.rcParams["axes.labelsize"]   = 17        #
  matplotlib.rcParams["xtick.labelsize"]  = 13        #
  matplotlib.rcParams["ytick.labelsize"]  = 13        #
  matplotlib.rcParams['image.cmap'] = 'plasma'        # Colormap default
  matplotlib.rcParams['axes.formatter.useoffset'] = False # another one I've forgotten the purpose of


def transect(data, z, lat, lon, start, end, npoints=100, 
             topog=None, latt=None, lont=None, ztop=4000,
             title="", ax=None, 
             cmap=None, norm=None, cbarform=None,
             contours=None,lines=None):
    '''
    Draw cross section
        data is 3d
        z (3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    # Potential temperature
    slicedata  = utils.cross_section(data,lat,lon,start,end,npoints=npoints)
    
    # Pull out cross section of topography and height
    if topog is not None:
        slicetopog = utils.cross_section(topog,latt,lont,start,end,npoints=npoints)
    
    slicez = utils.cross_section(z,lat,lon,start,end,npoints=npoints)
    #xticks,xlabels = utils.cross_section_ticks_labels(start,end)
    xaxis=np.linspace(0,1,npoints)
    slicex=np.tile(xaxis,(len(z),1))
    
    if ax is not None:
        plt.sca(ax)
    # Note that contourf can work with non-plaid coordinate grids provided both are 2-d
    # Contour inputs: xaxis, yaxis, data, colour gradient 
    if contours is None:
        plt.contourf(slicex,slicez,slicedata,cmap=cmap,norm=norm)
    else:
        plt.contourf(slicex,slicez,slicedata,contours,cmap=cmap,norm=norm)
    
    if cbarform is not None:
        plt.colorbar(format=cbarform)
    
    # Add contour lines
    if lines is not None:
        plt.contour(slicex,slicez,slicedata,lines,colors='k')            
    
    # make sure land is obvious
    if topog is not None:
        plt.fill_between(xaxis,slicetopog,interpolate=True,facecolor='black')
    
    if ztop != None:
        plt.ylim(0,ztop)
    
    plt.xticks([])
    plt.xlabel('')
    plt.title(title)

def transect_s(s, z, lat, lon, start, end, npoints=100, 
               topog=None, latt=None, lont=None, ztop=4000,
               title="Wind speed (m/s)", ax=None, 
               cmap=plt.cm.get_cmap('YlGnBu'), norm=None, cbarform=None,
               contours=np.arange(0,22,2),lines=np.arange(0,22,2)):
    '''
    Draw wind speed cross section
        s is 3d wind speed
        z(3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    
    # wind speed
    s[np.isnan(s)] = -5000 # There is one row or column of s that is np.NaN, one of the edges I think
    
    # call transect using some defaults for potential temperature
    transect(s,z,lat,lon,start,end,npoints=npoints,
             topog=topog, latt=latt, lont=lont, ztop=ztop,
             title=title, ax=ax, 
             cmap=cmap, norm=norm, cbarform=cbarform,
             contours=contours,lines=lines)

def transect_theta(theta, z, lat, lon, start, end, npoints=100, 
                   topog=None, latt=None, lont=None, ztop=4000,
                   title="$T_{\\theta}$ (K)", ax=None, 
                   cmap=plt.cm.get_cmap('YlOrRd'), norm=None, cbarform=None,
                   contours=np.arange(280,320,2),lines=np.arange(280,320,2)):
    '''
    Draw theta cross section
        theta is 3d potential temperature
        z(3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    # call transect using some defaults for potential temperature
    transect(theta,z,lat,lon,start,end,npoints=npoints,
             topog=topog, latt=latt, lont=lont, ztop=ztop,
             title=title, ax=ax, 
             cmap=cmap, norm=norm, cbarform=cbarform,
             contours=contours,lines=lines)

def transect_w(w, z, lat, lon, start, end, npoints=100, 
               topog=None, latt=None, lont=None, ztop=4000,
               title="Vertical motion (m/s)", ax=None, 
               cmap=plt.cm.get_cmap('PiYG'), norm=col.SymLogNorm(0.25), cbarform=tick.ScalarFormatter(),
               contours=np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0])),
               lines=np.array([0])):
    '''
    Draw theta cross section
        w is 3d vertical motion
        z(3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    if cmap is not None:
        cmap.set_over('k')
    # call transect using some defaults for vertical velocity w
    transect(w, z,lat,lon,start,end,npoints=npoints,
             topog=topog, latt=latt, lont=lont, ztop=ztop,
             title=title, ax=ax, 
             cmap=cmap, norm=norm, cbarform=cbarform,
             contours=contours,lines=lines)

def transect_qc(qc, z, lat, lon, start, end, npoints=100, 
               topog=None, latt=None, lont=None, ztop=4000,
               title="Water and Ice (kg/kg air)", ax=None, 
               cmap=plt.cm.get_cmap('YlGnBu'), norm=None, cbarform=None,
               contours=np.arange(0.0,2.0,0.25),
               lines=np.array([0.1])):
    '''
    Draw theta cross section
        w is 3d vertical motion
        z(3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    if cmap is not None:
        cmap.set_over('k')
    # call transect using some defaults for vertical velocity w
    transect(qc, z,lat,lon,start,end,npoints=npoints,
             topog=topog, latt=latt, lont=lont, ztop=ztop,
             title=title, ax=ax, 
             cmap=cmap, norm=norm, cbarform=cbarform,
             contours=contours,lines=lines)

def map_add_locations(namelist, proj=None,
                      marker='o', color='r', markersize=None, 
                      dx=.025,dy=.015):
    '''
    input is list of names to be added to a map using the lat lons in _latlons_
    '''
    for name in namelist:
        y,x=_latlons_[name]
        # Add marker
        plt.plot(x,y,  color=color, linewidth=0, 
                 marker=marker, markersize=None, transform=proj)
        # add name
        plt.text(x+dx, y+dy, name,
             horizontalalignment='right',
             transform=proj)
    

def map_google(extent,zoom=10,fig=None,subplotxyn=None,draw_gridlines=True,gridlines=None):
    '''
    Draw a map with google image tiles over given EWSN extent
    if this will be part of a larger figure, subplotxyn=[3,3,4] will put it there
    to define where gridlines go, put gridlines = [lats,lons] 
    return figure, axis, projection
    '''
    # Request map from google
    request = cimgt.GoogleTiles()
    gproj=request.crs
    # Use projection to set up plot
    if fig is None:
        fig = plt.figure()
    
    if subplotxyn is None:
        ax = fig.add_subplot(1, 1, 1, projection=gproj)#stamen_terrain.crs)
    else:
        nrows,ncols,n= subplotxyn
        ax = fig.add_subplot(nrows,ncols,n, projection=gproj)

    if draw_gridlines:
        print("drawing grid")
        gl = ax.gridlines(linewidth=1, color='black', alpha=0.5, linestyle='--', draw_labels=True)
        gl.xlabels_top = gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        if gridlines is not None:
            yrange,xrange=gridlines
            gl.xlocator = matplotlib.ticker.FixedLocator(xrange)
            gl.ylocator = matplotlib.ticker.FixedLocator(yrange)

    # Where are we looking
    ax.set_extent(extent)
    # default interpolation ruins location names
    ax.add_image(request, zoom, interpolation='spline36') 
    return fig, ax, gproj

def map_topography(extent, topog,lat,lon,title="Topography",cmap='copper'):
    '''
    Show topography map matching extents
    '''
    
    plt.contourf(lon,lat,topog, cmap=cmap)
    # set x and y limits to match extent
    xlims = extent[0:2] # East to West
    ylims = extent[2:] # South to North
    plt.ylim(ylims); plt.xlim(xlims)
    plt.colorbar(label="m")
    plt.title(title)
    ## Turn off the tick values
    plt.xticks([]); plt.yticks([])

def utm_from_lon(lon):
    """
    utm_from_lon - UTM zone for a longitude

    Not right for some polar regions (Norway, Svalbard, Antartica)

    :param float lon: longitude
    :return: UTM zone number
    :rtype: int
    """
    return np.floor( ( lon + 180 ) / 6) + 1

def scale_bar(ax, proj, length, location=(0.5, 0.05), linewidth=3,
              units='km', m_per_unit=1000):
    """
    http://stackoverflow.com/a/35705477/1072212
    ax is the axes to draw the scalebar on.
    proj is the projection the axes are in
    location is center of the scalebar in axis coordinates ie. 0.5 is the middle of the plot
    length is the length of the scalebar in km.
    linewidth is the thickness of the scalebar.
    units is the name of the unit
    m_per_unit is the number of meters in a unit
    """
    # find lat/lon center to find best UTM zone
    x0, x1, y0, y1 = ax.get_extent(proj.as_geodetic())
    #x0, x1, y0, y1 = ax.get_extent(proj)
    # Projection in metres
    utm = cartopy.crs.UTM(utm_from_lon((x0+x1)/2))
    # Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(utm)
    # Turn the specified scalebar location into coordinates in metres
    sbcx, sbcy = x0 + (x1 - x0) * location[0], y0 + (y1 - y0) * location[1]
    # Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbcx - length * m_per_unit/2, sbcx + length * m_per_unit/2]
    # buffer for scalebar
    buffer = [patheffects.withStroke(linewidth=5, foreground="w")]
    # Plot the scalebar with buffer
    ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k',
        linewidth=linewidth, path_effects=buffer)
    # buffer for text
    buffer = [patheffects.withStroke(linewidth=3, foreground="w")]
    # Plot the scalebar label
    t0 = ax.text(sbcx, sbcy, str(length) + ' ' + units, transform=utm,
        horizontalalignment='center', verticalalignment='bottom',
        path_effects=buffer, zorder=2)
    left = x0+(x1-x0)*0.05
    # Plot the N arrow
    #t1 = ax.text(left, sbcy, u'\u25B2\nN', transform=utm,
    #    horizontalalignment='center', verticalalignment='bottom',
    #    path_effects=buffer, zorder=2)
    # Plot the scalebar without buffer, in case covered by text buffer
    ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k',
        linewidth=linewidth, zorder=3)
