# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 14:56:00 2019
    Script dealing with map figures for publication
@author: jgreensl
"""

# plotting 
from matplotlib import colors, ticker, patches
import matplotlib.pyplot as plt
# maths
import numpy as np

# mapping library
import cartopy

# geometry to add onto cartopy
import cartopy.feature
import shapely.geometry
import cartopy.crs
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.transforms import offset_copy

# LOCAL IMPORTS
from utilities import utils, fio, constants, plotting

### 
## GLOBALS
###
_sn_ = 'maps' # scriptname

# waroona nest centre and extents for run1
waroona_lat0, waroona_lon0 = -32.9, 116.1
waroona_res = [[.036, 384], [.01, 576], [.0028,576]] # resolution, nlats for each nest
waroona_extents = [[waroona_lon0-rx*(nx//2), waroona_lon0+rx*(nx//2), waroona_lat0-rx*(nx//2), waroona_lat0+rx*(nx//2)] for (rx,nx) in waroona_res]
# same for sirivan run1:
sirivan_lat0, sirivan_lon0 = 

# Plot defaults
plotting.init_plots()

def waroona_nest():
    """
    show waroona model run nested grids on google map, resolution annotated
    """
    #f, ax, gproj = plotting.map_google(plotting._extents_['waroona'], zoom=11)
    
    # Australia:
    aust = [107,140,-40,-10]
    # Request map from google
    request = cimgt.GoogleTiles()
    gproj=request.crs
    # Use projection to set up plot
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=gproj)
    
    # Where are we looking
    ax.set_extent(aust)
    
    # default interpolation ruins location names
    zoom = 4
    ax.add_image(request, zoom, interpolation='spline36') 
    
    ## Show grid resolution maybe with annotations
    ## Add box around zoomed in area
    xy=waroona_lon0,waroona_lat0
    for i in range(3):
        res,nres = waroona_res[i]
        width=res*nres
        botleft = xy[0]-width/2.0, xy[1]-width/2.0
        
        ax.add_patch(patches.Rectangle(xy=botleft,
                                       width=width, 
                                       height=width,
                                       #facecolor=None,
                                       fill=False,
                                       edgecolor='red',
                                       linewidth=2,
                                       #linestyle='-',
                                       alpha=0.9,
                                       transform=cartopy.crs.PlateCarree()
                                       ))
        
        ## add text?
        ax.annotate(['Nest 1','Nest 2','N3'][i], xy=botleft, 
                    xycoords=cartopy.crs.PlateCarree()._as_mpl_transform(ax), color='k',
                    ha='left', va='top')
    
        ax.annotate('Nest %d: %dx%d %0.1f km squares'%(i+1, nres, nres, [3.5, 1.0, 0.3][i]),
                    xy=[.055, .95 - .05*i],
                    xycoords='axes fraction',
                    color='k', fontsize=13,
                    ha='left', va='top')

    ## Add zoom of nest 3?
    fio.save_fig('waroona_run1',_sn_,'nested_grid.png',plt)

def sirivan_nest():
    """
    show sirivan model run nested grids on google map, resolution annotated
    """
    #f, ax, gproj = plotting.map_google(plotting._extents_['waroona'], zoom=11)
    
    # East Australia:
    aust = [117,155,-40,-10]
    # Request map from google
    request = cimgt.GoogleTiles()
    gproj=request.crs
    # Use projection to set up plot
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=gproj)
    
    # Where are we looking
    ax.set_extent(aust)
    
    # default interpolation ruins location names
    zoom = 4
    ax.add_image(request, zoom, interpolation='spline36') 
    
    ## Show grid resolution maybe with annotations
    ## Add box around zoomed in area
    xy=plotting._latlons_['nest_centre'][::-1]
    for i in range(3):
        res,nres = plotting.__nest_res__[i]
        width=res*nres
        botleft = xy[0]-width/2.0, xy[1]-width/2
        
        ax.add_patch(patches.Rectangle(xy=botleft,
                                       width=width, 
                                       height=width,
                                       #facecolor=None,
                                       fill=False,
                                       edgecolor='red',
                                       linewidth=2,
                                       #linestyle='-',
                                       alpha=0.9,
                                       transform=cartopy.crs.PlateCarree()
                                       ))
        
        ## add text?
        ax.annotate(['Nest 1','Nest 2','N3'][i], xy=botleft, 
                    xycoords=cartopy.crs.PlateCarree()._as_mpl_transform(ax), color='k',
                    ha='left', va='top')
    
        ax.annotate('Nest %d: %dx%d %0.1f km squares'%(i+1, nres, nres, [3.5, 1.0, 0.3][i]),
                    xy=[.055, .95 - .05*i],
                    xycoords='axes fraction',
                    color='k', fontsize=13,
                    ha='left', va='top')

    ## Add zoom of nest 3?
    fio.save_fig('sirivan_run1',_sn_,'nested_grid.png',plt)