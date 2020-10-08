# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 11:16:28 2020

@author: Jesse
"""

# read tiff file
from osgeo import gdal

#
#import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patheffects

# Mapping stuff
#import cartopy
import cartopy.crs as ccrs
#from cartopy.io import shapereader
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#import cartopy.io.img_tiles as cimgt
#from owslib.wmts import WebMapTileService

def map_add_grid(ax, **gridargs):
    """
    add gridlines to some plot axis
    ARGUMENTS:
        ax: plt axis instance
        arguments for ax.gridlines
            defaults are draw_labels=True, color='gray', alpha=0.35, ls='--', lw=2
    """
    if 'draw_labels' not in gridargs:
        gridargs['draw_labels']=True
    if 'linewidth' not in gridargs:
        gridargs['linewidth']=2
    if 'color' not in gridargs:
        gridargs['color']='gray'
    if 'alpha' not in gridargs:
        gridargs['alpha']=0.35
    if 'linestyle' not in gridargs:
        gridargs['linestyle']='--'
    
    gl = ax.grid(gridargs)
    #gl.xlabels_top = False
    #gl.ylabels_left = False
    #gl.xlines = False
    
    return gl

def map_add_location(text, loc,
                     proj=None,
                     marker='o', color='grey', markersize=None, 
                     textcolor='k',
                     dx=.025,dy=.015):
    '''
    ARGS:
        text: label for marker
        loc: [lat,lon] of marker
        dx,dy: offset for textual label
    '''
    
    y,x=loc
    
    plotkwargs={"color":color,
                "linewidth":0,
                "marker":marker,
                "markersize":None,
                }
    
    textkwargs={"color":textcolor,
                "horizontalalignment":"right",
                }
    if proj is not None:
        plotkwargs["transform"]=proj
        textkwargs["transform"]=proj
    
    # Add marker and text
    plt.plot(x,y,  **plotkwargs)
    plt.text(x+dx, y+dy, text, **textkwargs)

def map_add_nice_text(ax, latlons, texts=None, markers=None, 
                      fontsizes=12, fontcolors='wheat', 
                      markercolors='grey', markersizes=None,
                      outlinecolors='k', transform=None):
    '''
    ARGUMENTS:
        ax: plot axis
        latlons: iterable of (lat, lon) pairs
        texts: iterable of strings, optional
        markers: iterable of characters, optional
        transform: if using geoaxes instance (non platecarree map) then set this to ccrs.Geodetic() or PlateCarree()?
    '''
    ## Adding to a map using latlons can use the geodetic transform
    #geodetic_CRS = ccrs.Geodetic()
    transformargs = {}
    if transform is not None:
        if transform == True:
            transform=ccrs.Geodetic()
        transformargs['transform']=transform
    
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
        ax.plot(lon, lat, color=mcolor, linewidth=0, 
                marker=marker, markersize=msize,
                path_effects=marker_effects,
                **transformargs)
        
        if len(text)>0:
            # Add text to map
            txt = ax.text(lon, lat, text, fontsize=fsize, color=fcolor,
                          **transformargs)
            # Add background (outline)
            txt.set_path_effects(text_effects)    

def map_tiff_qgis(fname='waroona_osm.tiff', extent=None, show_grid=False,
                  fig=None, subplot_row_col_n=[1,1,1], subplot_axes=None,
                  aspect='auto'):
    """
    fname needs to be on EPSG 4326 (platecarree)
    
    ARGUMENTS:
        fname: "sirivan_map_linescan.tiff" (Searches in data/QGIS/)
        extent: [left,right,bot,top] (lats and lons)
            where is projection zoomed in to
        show_grid: {True|False}
            show lats and lons around map
        fig: matplotlib pyplot figure (default=[1,1,1])
            can add map to an already created figure object
        subplot_row_col_n: [row,col,n] (optional)
            subplot axes to create and add map to
        subplot_extent: [x0,y0,width,height]
            where to put axes inside figure (can use instead of subplot_row_col_n)
        aspect: {'equal' | 'auto'} # default is 'auto', equal prevents pixel stretch
    """
    gdal.UseExceptions()
    
    #path_to_tiff = "data/QGIS/"+fname
    path_to_tiff = fname
    
    if subplot_row_col_n is None:
        subplot_row_col_n = [1,1,1]
    
    # gdal dataset
    # this prints a warning statement unless the .tiff has been transformed with tiff2rgba
    ds = gdal.Open(path_to_tiff)

    # RGBa image read in as a numpy array
    img = plt.imread(path_to_tiff)
    
    # projection defined in QGIS
    # 4326 is not actually 'projected' - except it is PlateCarree
    # projection=ccrs.PlateCarree() if EPSG == 4326 else ccrs.epsg(str(EPSG))
    
    # geotransform for tiff coords
    # tells us the image bounding coordinates
    gt = ds.GetGeoTransform()
    imgextent = (gt[0], gt[0] + ds.RasterXSize * gt[1],
                 gt[3] + ds.RasterYSize * gt[5], gt[3])
    
    ## set up figure and axes if needed
    if fig is None:
        fig = plt.figure(figsize=(13, 9))
    if subplot_axes is not None:
        ax = fig.add_axes(subplot_axes, frameon=False)
    else:
        nrows,ncols,n = subplot_row_col_n
        ax = fig.add_subplot(nrows,ncols,n)
    
    # Display tiff and cut to extent
    ax.imshow(img,
              extent=imgextent, 
              origin='upper',
              aspect=aspect)
    if extent is not None:
        #ax.set_extent(extent)
        ax.set_xlim([extent[0],extent[1]]) # E-W
        ax.set_ylim([extent[2],extent[3]]) # S-N
            
    if show_grid:
        map_add_grid(ax)
    
    return fig, ax


fig,ax = map_tiff_qgis("waroona_osm.tiff",)
# 
map_add_nice_text(ax,[[-32.84, 115.93],],texts=["WAROONA",])

plt.show()