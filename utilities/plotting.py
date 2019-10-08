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
# Extents: EWSN
_extents_               = {}
_latlons_               = {}

# Waroona locations
_extents_['waroona']    = [115.775,116.2, -33.05,-32.7] # local
_extents_['waroonas']   = [112,120,-34.5,-31] # synoptic
_extents_['waroonaz']    = [115.92, 116.19, -32.92,-32.83] # zoom in on fire
## Nests centre: -32.9, 116.1
## Nests resolution: 0.036 384x384, 0.01, 0.0028 (~300m)


_latlons_['waroona']    = -32.84, 115.93  # suburb centre: -32.8430, 115.8526
_latlons_['yarloop']    = -32.96, 115.90  # suburb centre: -32.9534, 115.9124
_latlons_['wagerup']    = -32.92, 115.91  # AWS at wagerup, 40 m asl
_latlons_['perth']      = -31.9505, 115.8605
_latlons_['fire_waroona'] = -32.89, 116.17
_latlons_['fire_waroona_upwind'] = -32.89 -0.004, 116.17+0.009 # ~ 1km from fire

_latlons_['nest_centre'] = -32.9, 116.1
__nest_res__ = [[.036, 384], [.01, 576], [.0028,576]] # resolution, nlats for each nest
for i in range(3):
    lat,lon = _latlons_['nest_centre']
    rx, nx = __nest_res__[i]
    hx=nx//2
    _extents_['waroona_nest%d'%(i+1)] = [lon-hx*rx, lon+hx*rx, lat-hx*rx, lat+hx*rx]

# two PyroCB
_latlons_['pyrocb_waroona1'] = -32.87,116.1 # ~4pm first day
_latlons_['pyrocb_waroona2'] = 0,0 # 1100-1400 second day

# Sir Ivan locations
_extents_['sirivan']    = [149.2, 150.4, -32.4, -31.6]
_extents_['sirivan_linescan']   = [149.48, 150.04, -32.18, -31.85]
_extents_['sirivanz']   = [149.4, 150.15, -32.2, -31.8]
_extents_['sirivans']   = [145,154, -34, -30]
_extents_['sirivan_linescans'] = [168,221, ]
_latlons_['dunedoo']    = -31.99, 149.53
_latlons_['uarbry']      = -32.047280, 149.71
_latlons_['cassillis']      = -32.01, 150.0
_latlons_['fire_sirivan'] = -32.45, 149.51
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
_transects_['pyrocb_waroona'] = [-32.89 , 115.9 ], [-32.88   , __x1__]
_transects_['pyrocbx1_waroona'] = [-32.95, 116.09], [-32.78   , 116.16]
_transects_['pyrocbx2_waroona'] = [-32.95, 116.15], [-32.78   , 116.09]
_transects_['emberstorm1'] = [-32.86, __x0__+.09], [-32.88, __x0__+.2] # emberstorm 

# looking at sir ivan
#_extents_['sirivan']    = [149.2, 150.4, -32.4, -31.6]
__si0__, __si1__ = 149.4, 149.9
_transects_['sirivan1'] = [-32.  , __si0__], [-32.3 , __si1__]
_transects_['sirivan2'] = [-32.05, __si0__], [-32.35, __si1__] # like sirivan1 but lower
_transects_['sirivan3'] = [-31.95, __si0__], [-32.25, __si1__]   # '' but higher
_transects_['sirivan4'] = [-32.25, __si0__], [-31.95, __si1__] # low lat to high lat crosssec
_transects_['sirivan5'] = [-32.35, __si0__], [-32.05, __si1__] 
_transects_['sirivan6'] = [-32.3 , __si0__], [-32   , __si1__] 




def init_plots():
  matplotlib.rcParams['font.size'] = 15.0
  matplotlib.rcParams["text.usetex"]      = False     # I forget what this is for, maybe allows latex labels?
  matplotlib.rcParams["legend.numpoints"] = 1         # one point for marker legends
  matplotlib.rcParams["figure.figsize"]   = (9, 7)    # Default figure size
  matplotlib.rcParams["axes.titlesize"]   = 19        # title font size
  matplotlib.rcParams["figure.titlesize"] = 21        # figure suptitle size
  matplotlib.rcParams["axes.labelsize"]   = 15        #
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
    if extentname == 'waroona':
        map_add_locations(['waroona','yarloop'], 
                          text=[['Waroona', 'Yarloop'],['','']][hide_text], 
                          textcolor='k')
        # add fire ignition
        map_add_locations(['fire_waroona'],
                          text = [['Fire ignition'],['']][hide_text], 
                          color='r', marker='*', 
                          textcolor='k')
        # add pyroCB
    elif extentname == 'waroonaz':
        map_add_locations(['waroona'], 
                          text=[['Waroona'],['']][hide_text], 
                          textcolor='k')
        # add fire ignition
        map_add_locations(['fire_waroona'],
                          text = [['Fire ignition'],['']][hide_text], 
                          color='r', marker='*', 
                          textcolor='k')
    else:
        map_add_locations(['dunedoo','cassillis','uarbry'],
                          text=[['Dunedoo','Cassillis','Uarbry'],['','']][hide_text],
                          dx=[-.02,-.02,.05], dy =[-.015,-.015,-.03],
                          textcolor='k')
        # add fire ignition
        map_add_locations(['fire_sirivan'],
                          text = [['Fire ignition']['']][hide_text], dx=.05,
                          color='r', marker='*', 
                          textcolor='k')

def transect(data, z, lat, lon, start, end, npoints=100, 
             topog=None, latt=None, lont=None, ztop=4000,
             title="", ax=None, colorbar=True,
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
        plt.contourf(slicex,slicez,slicedata,cmap=cmap,norm=norm,extend='max')
    else:
        plt.contourf(slicex,slicez,slicedata,contours,cmap=cmap,norm=norm,extend='max')
    
    if colorbar:
        plt.colorbar(format=cbarform, pad=0.01) # pad is distance from axes
    
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

    return slicedata, slicex, slicez

def transect_s(s, z, lat, lon, start, end, npoints=100, 
               topog=None, latt=None, lont=None, ztop=4000,
               title="Wind speed (m/s)", ax=None, colorbar=True,
               cmap='YlGnBu', norm=None, cbarform=None,
               contours=np.arange(0,25,2.5),lines=np.arange(0,25,2.5)):
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
    return transect(s,z,lat,lon,start,end,npoints=npoints,
                    topog=topog, latt=latt, lont=lont, ztop=ztop,
                    title=title, ax=ax, colorbar=colorbar,
                    cmap=cmap, norm=norm, cbarform=cbarform,
                    contours=contours,lines=lines)

def transect_theta(theta, z, lat, lon, start, end, npoints=100, 
                   topog=None, latt=None, lont=None, ztop=4000,
                   title="$T_{\\theta}$ (K)", ax=None, colorbar=True,
                   cmap='YlOrRd', norm=None, cbarform=None,
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
    return transect(theta,z,lat,lon,start,end,npoints=npoints,
                    topog=topog, latt=latt, lont=lont, ztop=ztop,
                    title=title, ax=ax, colorbar=colorbar,
                    cmap=cmap, norm=norm, cbarform=cbarform,
                    contours=contours,lines=lines)

def transect_w(w, z, lat, lon, start, end, npoints=100, 
               topog=None, latt=None, lont=None, ztop=4000,
               title="Vertical motion (m/s)", ax=None, colorbar=True,
               cmap=_cmaps_['verticalvelocity'] , norm=col.SymLogNorm(0.25), 
               cbarform=tick.ScalarFormatter(),
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
    #if cmap is not None:
    #    if hasattr(cmap,'set_over'):
    #        cmap.set_over('k')
    # call transect using some defaults for vertical velocity w
    return transect(w, z,lat,lon,start,end,npoints=npoints,
                    topog=topog, latt=latt, lont=lont, ztop=ztop,
                    title=title, ax=ax, colorbar=colorbar,
                    cmap=cmap, norm=norm, cbarform=cbarform,
                    contours=contours,lines=lines)

def transect_qc(qc, z, lat, lon, start, end, npoints=100, 
               topog=None, latt=None, lont=None, ztop=4000,
               title="Water and ice (g/kg air)", ax=None, colorbar=True,
               cmap=_cmaps_['qc'] , norm=col.SymLogNorm(0.02), cbarform=tick.ScalarFormatter(),
               contours=np.arange(0.0,0.4,0.01),
               lines=np.array([constants.cloud_threshold])):
    '''
    Draw theta cross section
        qc is 3d vertical motion
        z(3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    
    # call transect using some defaults for vertical velocity w
    return transect(qc, z,lat,lon,start,end,npoints=npoints,
                    topog=topog, latt=latt, lont=lont, ztop=ztop,
                    title=title, ax=ax, colorbar=colorbar,
                    cmap=cmap, norm=norm, cbarform=cbarform,
                    contours=contours,lines=lines)

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
        


def map_google(extent,zoom=10,fig=None,subplot_row_col_n=None, 
               subplot_extent=None, draw_gridlines=True,gridlines=None):
    '''
    Draw a map with google image tiles over given EWSN extent
    if this will be part of a larger figure, subplot_row_col_n=[3,3,4] will put it there
    to define where gridlines go, put gridlines = [lats,lons] 
    return figure, axis, projection
    '''
    # Request map from google
    request = cimgt.GoogleTiles()
    gproj=request.crs
    # Use projection to set up plot
    if fig is None:
        fig = plt.figure()
    
    if subplot_extent is not None:
        ax = fig.add_axes(subplot_extent, projection=gproj)
    elif subplot_row_col_n is not None:
        nrows,ncols,n= subplot_row_col_n
        ax = fig.add_subplot(nrows,ncols,n, projection=gproj)
    else:
        ax = fig.add_subplot(1,1,1, projection=gproj)

    if draw_gridlines:
        map_draw_gridlines(ax, gridlines=gridlines)

    # Where are we looking
    ax.set_extent(extent)
    # default interpolation ruins location names
    ax.add_image(request, zoom, interpolation='spline36') 
    return fig, ax, gproj

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
                 cbar=True, cbarform=None):
    '''
    Show topography map matching extents
    '''
    
    cs = plt.contourf(lon,lat,data, levels=clevs, cmap=cmap,norm=norm)
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

def map_satellite(extent = _extents_['waroona'], 
                  fig=None, subplot_row_col_n=None,
                  subplot_extent=None,
                  show_name=True, name_size=10):
    '''
    use NASA GIBS: Global Imagery Browse Services, to get high res satellite image: 
        https://wiki.earthdata.nasa.gov/display/GIBS/GIBS+Available+Imagery+Products#expand-SurfaceReflectance16Products
    
    '''
    URL = 'http://gibs.earthdata.nasa.gov/wmts/epsg4326/best/wmts.cgi'
    wmts = WebMapTileService(URL)
    
    layer = 'Landsat_WELD_CorrectedReflectance_TrueColor_Global_Annual'
    
    date_str = '2010-01-06' # landsat not available after 2010
    
    ## Plot setup
    # plot projections (for coord ref transforms)
    plot_CRS = ccrs.Mercator()
    geodetic_CRS = ccrs.Geodetic()
    
    # where are we looking?
    lon0, lon1, lat0, lat1 = extent
    # transform to map corners
    x0, y0 = plot_CRS.transform_point(lon0, lat0, geodetic_CRS)
    x1, y1 = plot_CRS.transform_point(lon1, lat1, geodetic_CRS)
    
    if fig is None:
        # keep aspect ratio based on lat lon corners
        ysize = 8
        xsize = ysize * (x1 - x0) / (y1 - y0)
        fig = plt.figure(figsize=(xsize, ysize), dpi=100)

    # create plot axes
    if subplot_extent is not None:
        ax = fig.add_axes(subplot_extent, projection=plot_CRS)
    elif subplot_row_col_n is not None:
        nrows,ncols,n= subplot_row_col_n
        ax = fig.add_subplot(nrows,ncols,n, projection=plot_CRS)
    else:
        ax = fig.add_subplot(1,1,1, projection=plot_CRS)
        
    ax.set_xlim((x0, x1))
    ax.set_ylim((y0, y1))
    
    ## add map tile from web service
    ax.add_wmts(wmts, layer, wmts_kwargs={'time': date_str})
    
    if show_name:
        # add layer name
        lat_bl = lat0 + 0.02*(lat1-lat0)
        lon_bl = lon0 + 0.05*(lon1-lon0)
        txt = ax.text(lon_bl, lat_bl, wmts[layer].title, 
                      fontsize=name_size, color='wheat',
                      transform=geodetic_CRS)
        txt.set_path_effects([patheffects.withStroke(linewidth=5,
                                                     foreground='black')])
    
    return fig, ax, plot_CRS, geodetic_CRS

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
        ax.plot(lon, lat, color=mcolor, linewidth=0, marker=marker, markersize=None,
                transform=geodetic_CRS,
                path_effects=marker_effects)
        
        if len(text)>0:
            # Add text to map
            txt = ax.text(lon, lat, text, fontsize=fsize, color=fcolor,
                          transform=geodetic_CRS)
            # Add background (outline)
            txt.set_path_effects(text_effects)
    
def map_topography(extent, topog,lat,lon,title="Topography"):
    '''
    Show topography map matching extents
    '''
    # push blue water part of scale a bit lower
    clevs= np.linspace(-150,550,50,endpoint=True)
    cmaptr=plt.cm.get_cmap("terrain")
    return map_contourf(extent, topog,lat,lon,title=title,clevs=clevs,
                        cmap=cmaptr,clabel="m")
    

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
