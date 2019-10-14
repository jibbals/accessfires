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
__NESTS__ = {'waroona_run1':{'centre':[-32.9, 116.1],
                             'resolution':[.036,.01,.0028],
                             'nlats':[384,576,576],
                             'nlons':[384,576,576],
                             'wider_view':[107,140,-40,-10],
                             },
            'sirivan_run1':{'centre':[-32.001400, 149.798600],
                            'resolution':[.036,.01,.0028],
                             'nlats':[384,576,576],
                             'nlons':[384,576,576],
                             'wider_view':[125,162,-44,-11],
                             },
            }

# Plot defaults
plotting.init_plots()

def show_nests(model_run='waroona_run1', annotate_res=True, title=''):
    """
    show nested grids on google map, resolution annotated
    """
    #f, ax, gproj = plotting.map_google(plotting._extents_['waroona'], zoom=11)
    nest=__NESTS__[model_run]
    
    
    # look at nests within wider view:
    aust = nest['wider_view']
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
    xy=nest['centre'][::-1]
    for i in range(3):
        # our model grids are square, same number of lats and lons
        res,nres = nest['resolution'][i], nest['nlats'][i]
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
        
        if annotate_res:
            ax.annotate('Nest %d: %dx%d %0.1f km squares'%(i+1, nres, nres, [3.5, 1.0, 0.3][i]),
                        xy=[.055, .95 - .05*i],
                        xycoords='axes fraction',
                        color='k', fontsize=13,
                        ha='left', va='top')
    if title=='':
        locname = str.capitalize(model_run.split('_')[0])
        if locname == 'Sirivan':
            locname = 'Sir Ivan'
        title = "Fire-Access model run nests at %s"%locname
    plt.title(title)
    ## Add zoom of nest 3?
    fio.save_fig(model_run,_sn_,'nested_grid',plt)
    
if __name__=='__main__':
    
    for mr in ['sirivan_run1','waroona_run1']:
        show_nests(mr)