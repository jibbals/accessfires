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
from matplotlib.projections import register_projection
from matplotlib.collections import LineCollection

# Mapping stuff
import cartopy
import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt
from owslib.wmts import WebMapTileService
import warnings
# read tiff file
from osgeo import gdal

# local stuff
from utilities import utils, constants

###
### GLOBALS
###

## Standard colour maps for use between figures
_cmaps_ = {}
_cmaps_['verticalvelocity'] = 'PiYG_r'  # jeff likes this for vert velocity
_cmaps_['windspeed']        = 'YlGnBu' # jeff likes this for wind speeds
_cmaps_['qc']               = 'BuPu' # white through to purple for water+ice content in air
_cmaps_['topog']            = 'terrain'
_cmaps_['th']               = 'plasma' # potential temperature
# Extents: EWSN
_extents_               = {}
_latlons_               = {}

# Waroona locations
_extents_['waroona']    = [115.775,116.2, -33.05,-32.7] # local
_extents_['waroonas']   = [112,120,-34.5,-31] # synoptic
_extents_['waroonaz']    = [115.88, 116.19, -32.92,-32.83] # zoom in on fire
## Nests centre: -32.9, 116.1
## Nests resolution: 0.036 384x384, 0.01, 0.0028 (~300m)


_latlons_['waroona']    = -32.84, 115.93  # suburb centre: -32.8430, 115.8526
_latlons_['yarloop']    = -32.96, 115.90  # suburb centre: -32.9534, 115.9124
_latlons_['wagerup']    = -32.92, 115.91  # wagerup
_latlons_['AWS_wagerup']    = -32.92, 115.91  # AWS at wagerup, 40 m asl
_latlons_['perth']      = -31.9505, 115.8605
_latlons_['fire_waroona'] = -32.89, 116.17
_latlons_['fire_waroona_upwind'] = -32.89 -0.004, 116.17+0.009 # ~ 1km from fire

# two PyroCB
_latlons_['pyrocb_waroona1'] = -32.87,116.1 # ~4pm first day
_latlons_['pyrocb_waroona2'] = 0,0 # 1100-1400 second day

# Sir Ivan locations
_extents_['sirivan']    = [149.2, 150.4, -32.4, -31.6]
_extents_['sirivan_linescan']   = [149.48, 150.04, -32.18, -31.85]
_extents_['sirivanz']   = [149.4, 150.15, -32.2, -31.8]
_extents_['sirivans']   = [147,154, -34, -29] # synoptic

_latlons_['dunedoo']    = -31.99, 149.53
_latlons_['uarbry']      = -32.047280, 149.71
_latlons_['sirivan']    = _latlons_['uarbry'] # no idea where sir ivan is..
_latlons_['cassillis']      = -32.01, 150.0
_latlons_['fire_sirivan'] = -32.05, 149.59
_latlons_['fire_sirivan_upwind'] = -32.01, 149.5
# one pyrocb
_latlons_['pyrocb_sirivan'] = 0,0 # 0530UTC=XXXX local time
# The pycb over Sir Ivan was around 0530 UTC on Sunday 12 Feb (or a bit after).
# The exact location is a bit tricky, because the pycb would be downstream of 
# the actual fire front, but around the location of Uarbry 
# (which was burned over) is probably a good place to start. 
# You'll probably have to move the cross section around a bit to see what 
# lat/longs get the most interesting picture.


_transects_             = {} 
__x0__,__x1__ = 115.8, 116.19

_transects_['waroona1'] = [-32.79   , __x0__], [-32.92   , __x1__]
_transects_['waroona2'] = [-32.82   , __x0__], [-32.93   , __x1__]
_transects_['waroona3'] = [-32.86   , __x0__], [-32.88   , __x1__]
# again but start more southerly and end more northerly
_transects_['waroona4'] = [-32.92   , __x0__], [-32.82   , __x1__]
_transects_['waroona5'] = [-32.96   , __x0__], [-32.85   , __x1__] 
_transects_['waroona6'] = [-32.87   , __x0__], [-32.89   , __x1__]

# looking at sir ivan
#_extents_['sirivan']    = [149.2, 150.4, -32.4, -31.6]
_transects_['sirivan1'] = [-32.05, 149.4  ], [-32.0  , 150.3 ]
_transects_['sirivan2'] = [-32.0 , 149.4  ], [-31.95 , 150.3 ]
_transects_['sirivan3'] = [-32.1 , 149.4  ], [-32.15 , 150.3 ]
_transects_['sirivan4'] = [-31.8 , 149.45 ], [-32.15 , 150.2 ]
_transects_['sirivan5'] = [-31.95, 149.4  ], [-31.80 , 150.2 ]
_transects_['sirivan6'] = [-31.7 , 149.5  ], [-32.1  , 150.3 ]


def init_plots():
    matplotlib.rcParams['font.size'] = 14.0
    matplotlib.rcParams["text.usetex"]      = False     # I forget what this is for, maybe allows latex labels?
    matplotlib.rcParams["legend.numpoints"] = 1         # one point for marker legends
    matplotlib.rcParams["figure.figsize"]   = (9, 7)    # Default figure size
    matplotlib.rcParams["axes.titlesize"]   = 17        # title font size
    matplotlib.rcParams["figure.titlesize"] = 20        # figure suptitle size
    matplotlib.rcParams["axes.labelsize"]   = 14        #
    matplotlib.rcParams["xtick.labelsize"]  = 11        #
    matplotlib.rcParams["ytick.labelsize"]  = 11        #
    matplotlib.rcParams['image.cmap'] = 'plasma'        # Colormap default
    matplotlib.rcParams['axes.formatter.useoffset'] = False # another one I've forgotten the purpose of
    # rcParams["figure.dpi"] 
    #matplotlib.rcParams["figure.dpi"] = 400           # DEFAULT DPI for plot output
    # THIS MESSES UP THE PLOTTING


def map_add_locations_extent(extentname, hide_text=False):
    '''
    wrapper for map_add_locations that adds all the points for that extent
    '''
    
    locstrings = {'waroona':['waroona','yarloop'],
                  'waroonaz':['waroona'],
                  'sirivan':['dunedoo','cassillis','uarbry'],
                  'sirivans':['dunedoo','cassillis','uarbry']}
    dx=.025
    dxfire = .025
    dy=.015
    dyfire = .015
    if extentname in ['sirivan','sirivans']:
        dx=[.065,.02,.125]
        dy =[.02,.015,-.07]
        dxfire=.05
        dyfire=-.06

    # Where is fire ignition
    firename = "fire_"+extentname
    if extentname[-1] =="s":
        firename = "fire_"+extentname[:-1]

    locs = locstrings[extentname]
    text = [name.capitalize() for name in locs]
    if hide_text:
        text = ['']*len(locs)
    
    map_add_locations(locs, text=text, textcolor='k', dx=dx,dy=dy)
    # add fire ignition
    map_add_locations([firename], text=[['Ignition',''][hide_text]], 
                      color='r', marker='*', dx=dxfire, dy=dyfire, textcolor='k')
    # add weather stations
    if extentname=='waroona':
        map_add_locations(['wagerup'],[['AWS',''][hide_text]],
                          color='b', marker='*', dx=-.025,dy=.01)

def map_add_locations(namelist, text=None, proj=None,
                      marker='o', color='grey', markersize=None, 
                      textcolor='k',
                      dx=.025,dy=.015):
    '''
    input is list of names to be added to a map using the lat lons in _latlons_
    '''
    for i,name in enumerate(namelist):
        y,x=_latlons_[name]
        # maybe want special text
        if text is not None:
            name = text[i]
        # dx,dy can be scalar or list
        dxi,dyi = dx,dy
        if isinstance(dx, (list,tuple,np.ndarray)):
            dxi,dyi = dx[i], dy[i]

        # Add marker and text
        if proj is None:
            plt.plot(x,y,  color=color, linewidth=0, 
                     marker=marker, markersize=None)
            plt.text(x+dxi, y+dyi, name, color=textcolor,
                     horizontalalignment='right')
        else:
            plt.plot(x,y,  color=color, linewidth=0, 
                     marker=marker, markersize=None, transform=proj)
            plt.text(x+dxi, y+dyi, name, color=textcolor,
                     horizontalalignment='right',
                     transform=proj)


def map_draw_gridlines(ax, linewidth=1, color='black', alpha=0.5, 
                       linestyle='--', draw_labels=True,
                       gridlines=None):
    gl = ax.gridlines(linewidth=linewidth, color=color, alpha=alpha, 
                      linestyle=linestyle, draw_labels=draw_labels)
    if draw_labels:
        gl.xlabels_top = gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
    if gridlines is not None:
        yrange,xrange=gridlines
        gl.xlocator = matplotlib.ticker.FixedLocator(xrange)
        gl.ylocator = matplotlib.ticker.FixedLocator(yrange)

def map_contourf(extent, data, lat,lon, title="",
                 cmap=None, clabel="", clevs=None, norm=None, 
                 cbar=True, cbarform=None, **contourfargs):
    '''
    Show topography map matching extents
    '''
    
    cs = plt.contourf(lon,lat,data, levels=clevs, cmap=cmap,norm=norm, **contourfargs)
    cb = None
    
    # set x and y limits to match extent
    xlims = extent[0:2] # East to West
    ylims = extent[2:] # South to North
    plt.ylim(ylims); plt.xlim(xlims)
    if cbar:
        cb=plt.colorbar(label=clabel, format=cbarform, pad=0.01)
    plt.title(title)
    ## Turn off the tick values
    plt.xticks([]); plt.yticks([])
    return cs, cb

def map_fire(ff,lats,lons):
    """   """
    # only plot if there is fire
    if np.sum(ff<0) > 0:
        plt.contour(lons,lats,np.transpose(ff),np.array([0]), colors='red',
                    transform=ccrs.PlateCarree())


def map_tiff_qgis(file='sirivan.tiff', extent=None, show_grid=False,
                  locnames=None,
                  fig=None, subplot_row_col_n=[1,1,1], subplot_axes=None):
    """
    satellite image from ESRI, roads added from OSM using QGIS, saved as .tiff
    
    ARGUMENTS:
        locname: "sirivan_map_linescan.tiff" or whatever in qgis folder
        extent: [lon0, lon1, lat0, lat1]
            where is projection zoomed in to
        show_grid: {True|False}
            show lats and lons around map
        fig: matplotlib pyplot figure (default=[1,1,1])
            can add map to an already created figure object
        subplot_row_col_n: [row,col,n] (optional)
            subplot axes to create and add map to
        subplot_extent: [x0,y0,width,height]
            where to put axes inside figure (can use instead of subplot_row_col_n)
    """
    gdal.UseExceptions()
    
    path_to_tiff = "data/QGIS/"+file
    
    if subplot_row_col_n is None:
        subplot_row_col_n = [1,1,1]
    
    # gdal dataset
    # this prints a warning statement... pipe output and warning filter don't work
    ds = gdal.Open(path_to_tiff)

    # RGBa image read in as a numpy array
    img = plt.imread(path_to_tiff)
    
    # projection defined in QGIS
    projection = ccrs.epsg("3857")
    
    # geotransform for tiff coords
    # tells us the image bounding coordinates
    gt = ds.GetGeoTransform()
    imgextent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
                 gt[3] + ds.RasterYSize * gt[5], gt[3])
    
    if fig is None:
        fig = plt.figure(figsize=(13, 9))
    
    # if the subplotaxes are defined, use them, otherwise use subplot_row_col
    
    if subplot_axes is not None:
        ax = fig.add_axes(subplot_axes, frameon=False, projection=projection)
    else:
        nrows,ncols,n = subplot_row_col_n
        ax = fig.add_subplot(nrows,ncols,n,projection=projection)
    
    ax.imshow(img,
              extent=imgextent, 
              origin='upper')
    
    if extent is not None:
        ax.set_extent(extent)
        #crop to desired extent
        
        # if we have locations and defined extent, put them in 
        if locnames is not None:
            for locstr in locnames:
                loc=_latlons_[locstr]
                if (loc[0]<extent[2]) or (loc[0]>extent[3]) or (loc[1]<extent[0]) or (loc[1]>extent[1]):
                    locnames.remove(locstr)
        
    if show_grid:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=2, color='gray', alpha=0.35, linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_left = False
        #gl.xlines = False
    if locnames is not None:
        # split out fire into it's own thing
        fires = ['fire' in element for element in locnames]
        locations=[_latlons_[locstr] for locstr in locnames]
        LocationsNames = [str.capitalize(locstr) for locstr in locnames]
        markers=['o']*len(locations)
        markercolors=['grey']*len(locations)
        
        for i in range(len(locations)):
            if fires[i]:
                LocationsNames[i]=''
                markers[i] = '*'
                markercolors[i] = 'r'
        
        map_add_nice_text(ax, locations,
                          LocationsNames,
                          markers=markers,
                          markercolors=markercolors,
                          fontsizes=14)        
        
    return fig, ax, projection

def map_tiff_gsky(locname='waroona', fig=None, subplot_row_col_n=None,
             extent=None, show_grid=False, locnames=None):
    """
    satellite image from gsky downloaded into data folder, 
    used as background for figures
    
    ARGUMENTS:
        locname: { 'waroona' | 'sirivan' }
            which tiff are we using
        fig: matplotlib pyplot figure (default=[1,1,1])
            can add map to an already created figure object
        subplot_row_col_n: [row,col,n] (optional)
            subplot axes to create and add map to
        extent: [lon0, lon1, lat0, lat1]
            where is projection zoomed in to
        show_grid: {True|False}
            show lats and lons around map
        add_locations: {True|False}
            show nice text and markers for main points of interest
    """
    gdal.UseExceptions()
    
    path_to_tiff = "data/Waroona_Landsat_8_2015.tiff"
    
    if locname=='sirivan':
        path_to_tiff = "data/Sirivan_Landsat_8_2016.tiff"
    
    elif locname=='waroona_big':
        path_to_tiff = "data/waroona_big.tiff"
    
    
    if subplot_row_col_n is None:
        subplot_row_col_n = [1,1,1]
    
    # units are degrees latitude and longitude, so platecarree should work...
    projection = ccrs.PlateCarree()
    
    # gdal dataset
    ds = gdal.Open(path_to_tiff)
    # ndarray of 3 channel raster (?) [3, 565, 1024]
    data = ds.ReadAsArray()
    # np.sum(data<0) # 9 are less than 0
    data[data<0] = 0
    # 3 dim array of red green blue values (I'm unsure on what scale)
    R,G,B = data[0], data[1], data[2]
    # set them to a float array from 0 to 1.0
    scaled0 = np.array([R/np.max(R), G/np.max(G), B/np.max(B)])
    # image is too dark
    # gain=.3 and bias=.2 to increase contrast and brightness
    scaled = scaled0*1.3 + 0.1
    scaled[scaled>1] = 1.0
    # geotransform for tiff coords (?)
    gt = ds.GetGeoTransform()
    # extent = lon0, lon1, lat0, lat1
    fullextent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
                  gt[3] + ds.RasterYSize * gt[5], gt[3])
    
    ## Third: the figure
    if fig is None:
        fig = plt.figure(figsize=(13, 9))
    #subplot_kw = dict(projection=projection)
    
    # create plot axes
    nrows, ncols, n = subplot_row_col_n
    ax = fig.add_subplot(nrows, ncols, n, projection=projection)
    
    ax.imshow(scaled[:3, :, :].transpose((1, 2, 0)), 
              extent=fullextent,
              origin='upper')
    
    if extent is not None:
        ax.set_extent(extent)
        # Remove any locations outside of the extent
        if locnames is not None:
            for locstr in locnames:
                loc=_latlons_[locstr]
                if (loc[0]<extent[2]) or (loc[0]>extent[3]) or (loc[1]<extent[0]) or (loc[1]>extent[1]):
                    locnames.remove(locstr)
        
    if show_grid:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linewidth=2, color='gray', alpha=0.35, linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_left = False
        #gl.xlines = False
    if locnames is not None:
        # split out fire into it's own thing
        fires = ['fire' in element for element in locnames]
        locations=[_latlons_[locstr] for locstr in locnames]
        LocationsNames = [str.capitalize(locstr) for locstr in locnames]
        markers=['o']*len(locations)
        markercolors=['grey']*len(locations)
        
        for i in range(len(locations)):
            if fires[i]:
                LocationsNames[i]=''
                markers[i] = '*'
                markercolors[i] = 'r'
        
        map_add_nice_text(ax, locations,
                          LocationsNames,
                          markers=markers,
                          markercolors=markercolors,
                          fontsizes=14)        
        
    return fig, ax, projection


def map_quiver(u, v, lats, lons, nquivers=13, **quiver_kwargs):
    """
    wrapper for nice quiver overlay on a lat/lon map
    """
    # just want a fixed number of quivers
    nlats, nlons = u.shape
    latinds = np.round(np.linspace(0,nlats-1,nquivers)).astype(np.int)
    loninds = np.round(np.linspace(0,nlons-1,nquivers)).astype(np.int)
    u0 = u[latinds,:]
    u0 = u0[:,loninds]
    v0 = v[latinds,:]
    v0 = v0[:,loninds]
    # subset lats and lons
    if len(np.shape(lats))==1:
        sublats=lats[latinds]
        sublons=lons[loninds]
    # maybe we want to use 2D lats/lons
    elif len(np.shape(lats))==2:
        sublats=lats[latinds,:]
        sublats=sublats[:,loninds]
        sublons = lons[latinds,:]
        sublons = sublons[:,loninds]
    else:
        print("ERROR: Why are the lats 3 dimensional?")
        assert False
        
    ## Default quiver scale:
    if 'scale' not in quiver_kwargs:
        quiver_kwargs['scale'] = 85
        
    ## Default quiver pivot:
    if 'pivot' not in quiver_kwargs:
        quiver_kwargs['pivot'] = 'middle'
    
    plt.quiver(sublons, sublats, u0, v0, **quiver_kwargs)

def map_satellite(extent = _extents_['waroona'], 
                  fig=None, subplot_row_col_n=None,
                  subplot_extent=None,
                  show_name=True, name_size=10):
    '''
    Originally this used an API to pull some satellite imagery
    this proved unworkable through various bom and nci firewalls
    now this points to the tiff method
    '''
    locname='waroona'
    if extent[0] > 130:
        locname='sirivan'
        
    return map_tiff_gsky(locname, fig=fig, subplot_row_col_n=subplot_row_col_n,
                         extent=subplot_extent)

def map_add_nice_text(ax, latlons, texts=None, markers=None, 
                      fontsizes=12, fontcolors='wheat', 
                      markercolors='grey', markersizes=None,
                      outlinecolors='k'):
    '''
    ARGUMENTS:
        ax: plot axis
        latlons: iterable of (lat, lon) pairs
        texts: iterable of strings, optional
        markers: iterable of characters, optional
    '''
    ## Adding to a map using latlons can use the geodetic transform
    geodetic_CRS = ccrs.Geodetic()
    
    # Make everything iterable
    if texts is None:
        texts = ''*len(latlons)
    if markers is None:
        markers = 'o'*len(latlons)
    if isinstance(fontsizes, (int,float)):
        fontsizes = [fontsizes]*len(latlons)
    if (markersizes is None) or (isinstance(markersizes, (int,float))):
        markersizes = [markersizes]*len(latlons)
    if isinstance(fontcolors, str):
        fontcolors = [fontcolors]*len(latlons)
    if isinstance(markercolors, str):
        markercolors = [markercolors]*len(latlons)
    if isinstance(outlinecolors, str):
        outlinecolors = [outlinecolors]*len(latlons)

    
    
    # add points from arguments:
    for (lat, lon), text, marker, mcolor, msize, fcolor, fsize, outlinecolor in zip(latlons, texts, markers, markercolors, markersizes, fontcolors, fontsizes, outlinecolors):
        
        # outline for text/markers
        marker_effects=[patheffects.Stroke(linewidth=5, foreground=outlinecolor), patheffects.Normal()]
        text_effects = [patheffects.withStroke(linewidth=3, foreground=outlinecolor)]
        
        # Add the point to the map with patheffect
        ax.plot(lon, lat, color=mcolor, linewidth=0, marker=marker, markersize=msize,
                transform=geodetic_CRS,
                path_effects=marker_effects)
        
        if len(text)>0:
            # Add text to map
            txt = ax.text(lon, lat, text, fontsize=fsize, color=fcolor,
                          transform=geodetic_CRS)
            # Add background (outline)
            txt.set_path_effects(text_effects)
    
def map_topography(extent, topog,lat,lon,title="Topography", cbar=True):
    '''
    Show topography map matching extents
    '''
    # push blue water part of scale a bit lower
    clevs = np.linspace(-150,550,50,endpoint=True)
    # sir ivan fire is at higher altitudes
    if extent[0] > 140:
        clevs = np.linspace(100,800,50,endpoint=True)
    cmaptr=plt.cm.get_cmap("terrain")
    return map_contourf(extent, topog, lat, lon, 
                        title=title, clevs=clevs, cmap=cmaptr, 
                        clabel="m", cbar=cbar, cbarform=tick.ScalarFormatter())

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

def transect(data, z, lat, lon, start, end, npoints=100, 
             topog=None, ff=None, latt=None, lont=None, ztop=4000,
             title="", ax=None, colorbar=True,
             cbarform=None, contours=None,lines=None,
             **contourfargs):
    '''
    Draw cross section
        data is 3d
        z (3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    ## Default contourfargs
    if 'extend' not in contourfargs:
        contourfargs['extend'] = 'max'
    
    ## Check that z includes topography (within margin of 40 metres)
    if np.mean(z[0]+40)<np.mean(topog):
        print("ERROR:",np.mean(z[0]), np.min(z[0]), "(mean,lowest z) is lower than topog", np.mean(topog), np.min(topog))
        print("ERROR:", "Try adding topog to each level of z")
        assert False
        
    # Potential temperature
    slicedata  = utils.cross_section(data,lat,lon,start,end,npoints=npoints)
    
    # Pull out cross section of topography and height
    if latt is None:
        latt=lat
    if lont is None:
        lont=lon
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
        plt.contourf(slicex,slicez,slicedata,**contourfargs)
    else:
        plt.contourf(slicex,slicez,slicedata,contours,**contourfargs)
    
    if colorbar:
        plt.colorbar(format=cbarform, pad=0.01) # pad is distance from axes
    
    # Add contour lines
    if lines is not None:
        with warnings.catch_warnings():
            # ignore warning when there is no fire:
            warnings.simplefilter('ignore')
            plt.contour(slicex,slicez,slicedata,lines,colors='k')            
    
    # make sure land is obvious
    if topog is not None:
        plt.fill_between(xaxis,slicetopog,interpolate=True,facecolor='darkgrey')
    
    # put it red where the fire is burnt
    if ff is not None:
        ffslice = utils.cross_section(ff, lat, lon, start, end, 
                                      npoints=npoints)
        burnt = np.copy(slicetopog)
        interp=False
        # 0.005 is close enough?
        if not np.all(ffslice>=0.005):
            burnt[ffslice>=0.005] = np.NaN
            plt.fill_between(xaxis, burnt, interpolate=interp, 
                             facecolor='yellow', zorder=1)
        # 0.003 is pretty near the fire front
        if not np.all(ffslice>=0.003):
            burnt[ffslice>=0.003] = np.NaN
            plt.fill_between(xaxis, burnt, interpolate=interp, 
                             facecolor='orange', zorder=2)
        # 0 is the front
        if not np.all(ffslice>=0.001):
            burnt[ffslice>=0.001] = np.NaN
            plt.fill_between(xaxis, burnt, interpolate=interp, 
                             facecolor='red', zorder=3)
        # -ve is behind the front
        if not np.all(ffslice>=0):
            burnt[ffslice>0]=np.NaN
            plt.fill_between(xaxis, burnt, interpolate=interp, 
                             facecolor='saddlebrown', zorder=4)
    
    if ztop != None:
        plt.ylim(0,ztop)
    
    plt.xticks([])
    plt.xlabel('')
    plt.title(title)

    return slicedata, slicex, slicez

def transect_s(s, z, lat, lon, start, end, npoints=100, 
               topog=None, ff=None, latt=None, lont=None, ztop=4000,
               title="Wind speed (m/s)", ax=None, colorbar=True,
               cbarform=None, contours=np.arange(0,25,2.5),
               lines=np.arange(0,25,2.5), **contourfargs):
    '''
    Draw wind speed cross section
        s is 3d wind speed
        z(3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    ## default cmap
    if 'cmap' not in contourfargs:
        contourfargs['cmap']=_cmaps_['windspeed']
    
    # wind speed
    s[np.isnan(s)] = -5000 # There is one row or column of s that is np.NaN, one of the edges I think
    
    # call transect using some defaults for potential temperature
    return transect(s,z,lat,lon,start,end,npoints=npoints,
                    topog=topog, ff=ff, latt=latt, lont=lont, ztop=ztop,
                    title=title, ax=ax, colorbar=colorbar,
                    cbarform=cbarform, contours=contours,lines=lines,
                    **contourfargs)

def transect_theta(theta, z, lat, lon, start, end, npoints=100, 
                   topog=None, ff=None, latt=None, lont=None, ztop=4000,
                   title="$T_{\\theta}$ (K)", ax=None, colorbar=True,
                   cbarform=tick.ScalarFormatter(),
                   contours = np.arange(280,350,1),
                   lines = np.union1d(np.arange(280,301,2), np.arange(310,351,10)),
                   **contourfargs):
    '''
    Draw theta cross section
        theta is 3d potential temperature
        z(3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    if 'cmap' not in contourfargs:
        contourfargs['cmap'] = _cmaps_['th']
    if 'norm' not in contourfargs:
        contourfargs['norm'] = col.SymLogNorm(300)
        
    # call transect using some defaults for potential temperature
    return transect(theta,z,lat,lon,start,end,npoints=npoints,
                    topog=topog, ff=ff, latt=latt, lont=lont, ztop=ztop,
                    title=title, ax=ax, colorbar=colorbar,
                    cbarform=cbarform, contours=contours, lines=lines,
                    **contourfargs)

def transect_w(w, z, lat, lon, start, end, npoints=100, 
               topog=None, ff=None, latt=None, lont=None, ztop=4000,
               title="Vertical motion (m/s)", ax=None, colorbar=True, 
               cbarform=tick.ScalarFormatter(),
               contours=np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0])),
               lines=np.array([0]),
               **contourfargs):
    '''
    Draw theta cross section
        w is 3d vertical motion
        z(3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    if 'cmap' not in contourfargs:
        contourfargs['cmap'] = _cmaps_['verticalvelocity']
    if 'norm' not in contourfargs:
        contourfargs['norm'] = col.SymLogNorm(0.25)
    
    # call transect using some defaults for vertical velocity w
    return transect(w, z,lat,lon,start,end,npoints=npoints,
                    topog=topog, ff=ff, latt=latt, lont=lont, ztop=ztop,
                    title=title, ax=ax, colorbar=colorbar,
                    cbarform=cbarform, contours=contours, lines=lines,
                    **contourfargs)

def transect_qc(qc, z, lat, lon, start, end, npoints=100, 
               topog=None, ff=None, latt=None, lont=None, ztop=4000,
               title="Water and ice (g/kg air)", ax=None, colorbar=True,
               cbarform=tick.ScalarFormatter(), contours=np.arange(0.0,0.4,0.01),
               lines=np.array([constants.cloud_threshold]),
               **contourfargs):
    '''
    Draw theta cross section
        qc is 3d vertical motion
        z(3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    # defaults for contourfargs
    if 'cmap' not in contourfargs:
        contourfargs['cmap'] = _cmaps_['qc']
    if 'norm' not in contourfargs:
        contourfargs['norm'] = col.SymLogNorm(0.02)
    # call transect using some defaults for vertical velocity w
    return transect(qc, z,lat,lon,start,end,npoints=npoints,
                    topog=topog, ff=ff, latt=latt, lont=lont, ztop=ztop,
                    title=title, ax=ax, colorbar=colorbar,
                    cbarform=cbarform, contours=contours, lines=lines,
                    **contourfargs)

def transect_ticks_labels(start,end):
  '''
    return xticks and xlabels for a cross section
  '''
  lat1,lon1=start
  lat2,lon2=end
  # Set up a tuple of strings for the labels. Very crude!
  xticks = (0.0,0.5,1.0)
  fmt = '{:.1f}S {:.1f}E'
  xlabels = (fmt.format(-lat1,lon1),fmt.format(-0.5*(lat1+lat2),0.5*(lon1+lon2)),fmt.format(-lat2,lon2))
  return xticks,xlabels
 

def utm_from_lon(lon):
    """
    utm_from_lon - UTM zone for a longitude

    Not right for some polar regions (Norway, Svalbard, Antartica)

    :param float lon: longitude
    :return: UTM zone number
    :rtype: int
    """
    return np.floor( ( lon + 180 ) / 6) + 1


#### NOT WORKING/FINISHED ####
# def transect_wall(model_run = 'waroona_run3', 
#                 hour=20,
#                 extent=[115.88, 116.0, -32.86,-32.83],
#                 HSkip=None,
#                 send_to_browser=True):
#     """
#         look at waroona emberstorm bits
#             surface wind speeds
#             heat flux
#             pot temp up to ~ 600m
#         TODO: 
#             make transect show higher altitude, but cap whole figure at some height
#             add surface wind direction/speed somehow
#             vertical motion instead of pot temp on transect wall
#             move transect wall to middle of figure, make it .8 opaque
            
#     """
#     top_height=1000
#     hours=fio.model_outputs[model_run]['filedates']
#     dtime=hours[hour]
#     # [lat,lon], [lat1,lon1] of extent transect through middle
#     transect = [[(extent[2]+extent[3])/2,extent[0]],[(extent[2]+extent[3])/2,extent[1]]]
    
#     cubes = fio.read_model_run(model_run, 
#                                fdtime=dtime, 
#                                extent=extent, 
#                                add_theta=True, 
#                                add_topog=True, 
#                                add_winds=True,
#                                add_z=True,
#                                HSkip=HSkip)
    
#     th, qc, z_th = cubes.extract(['potential_temperature','qc','z_th'])
#     # datetimes in hour output
#     cubetimes = utils.dates_from_iris(th)
    
#     ff, sh = fio.read_fire(model_run,
#                            dtimes=cubetimes,
#                            extent=extent,
#                            sensibleheat=True,
#                            HSkip=HSkip)
    
#     # Get the rest of the desired data
#     topog, = cubes.extract(['surface_altitude'])
#     d_topog = topog.data.data.T # convert to lon,lat
#     u,v = cubes.extract(['u','v'])
#     levh  = qc.coord('level_height').points
#     topind = np.sum(levh<top_height)
    
#     # these are level, lat, lon cubes
#     lat,lon = th.coord('latitude').points, th.coord('longitude').points
    
#     # dimensional mesh
#     X,Y,Z = np.meshgrid(lon,lat,levh) 
#     ## X Y Z are now [lat, lon, lev] for some reason
#     [X,Y,Z] = [ np.moveaxis(arr,0,1) for arr in [X,Y,Z]]
#     ## Now they are lon, lat, lev
#     ## Cut down to desired level
#     [X, Y, Z] = [ arr[:,:,:topind] for arr in [X,Y,Z]]
    
    
#     for ti, cubetime in enumerate(cubetimes):
#         surface_list=[]
        
#         # surface sensible heat flux [t,lon,lat]
#         d_sh = sh[ti].data.data
        
#         d_zth = z_th[:topind,:,:].data.data # [lev, lat, lon]
#         d_th = th[ti,:topind,:,:].data.data # [lev, lat, lon]
#         #d_th = threedee.cube_to_xyz(th[ti],ztopind=topind) # [lon,lat,lev]
        
#         # topography surface
#         topog_layer = go.Surface(
#             z=d_topog,
#             x=X[:,:,0],
#             y=Y[:,:,0],
#             colorscale=shcolor, # surface heat colormap 
#             cmin=shcmin,
#             cmax=shcmax,
#             #reversescale=True,
#             surfacecolor=d_sh, # colour by sensible heat
#             #opacityscale=[[0.0, 0], [100.0, .8], [shcmax, 1]], 
#             #hidesurface=True,
#             #showscale=False, # remove colour bar,
#             coloraxis='coloraxis', # first colour bar
#         )
#         surface_list.append(topog_layer)
        
#         ## Pull out transect to paint against the wall
#         #print("DEBUG:",d_th.shape, np.shape(lat),np.shape(lon),transect)
#         # transect of pot temp
#         sliceth  = utils.cross_section(d_th,lat,lon,transect[0],transect[1],npoints=len(lon))
#         # transect of topography
#         slicetopog = utils.cross_section(d_topog.T,lat,lon,transect[0],transect[1],npoints=len(lon))
#         # transect of height
#         slicez = utils.cross_section(d_zth,lat,lon,transect[0],transect[1],npoints=len(lon))
#         # transects are [lat, lon] - and should represent one latitude over multiple longitudes
#         #print(np.shape(slicez), np.shape(sliceth))
#         transect_layer = go.Surface(
#             z=slicez.T, # height
#             x=X[:,:,0], # longitudes
#             y=np.zeros(np.shape(slicez.T)) + np.max(Y), # highest latitude
#             colorscale='plasma', # pot temp colour map
#             cmin=380,
#             cmax=410,
#             #reversescale=True,
#             surfacecolor=sliceth.T, # colour by pot temp
#             #opacityscale=[[0.0, 0], [100.0, .8], [shcmax, 1]], 
#             #hidesurface=True,
#             #showscale=False, # remove colour bar,
#             coloraxis='coloraxis2', # second colour bar
#             colorbar = dict(xanchor='right'),
#         )
#         surface_list.append(transect_layer)
        
        
#         ### Paint figure
#         figname = None
#         if not send_to_browser:
#             #figname = cubetime.strftime('figures/threedee/test_%Y%m%d%H%M.png')
#             figname = fio.standard_fig_name(model_run,_sn_,cubetime,subdir='downslope')
        
#         layoutargs = dict(title=cubetime.strftime('%Y%m%d%H%M(UTC)'),
#                           font=dict(size=18, color="#111111"),)
                                    
#         create_figure(surface_list, filename=figname, **layoutargs)
