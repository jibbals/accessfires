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

##############################################
############### METHODS ####################
##############################################

def outline_waroona():
    """
    Plot tiff image of waroona, within southern WA.
        Also show zoomed in to fire area, and topography
    """
    
    extentname='waroonas'
    extent = plotting._extents_[extentname]
    inner = plotting._extents_['waroona']
    
    plotting.init_plots()
    
    cube_topog = fio.read_topog('waroona_old',inner)
    latt = cube_topog.coord('latitude').points
    lont = cube_topog.coord('longitude').points
    topog = cube_topog.data.data
    
    # Google map image tiles view of synoptic map
    fig, ax, proj = plotting.map_tiff(locname='waroona',
                                      extent=extent, #fig=fig,
                                      subplot_row_col_n=[2,1,1],
                                      show_grid=False, add_locations=True)
    #    fig,ax,proj=plotting.map_google(extent,
    #                                    zoom=6,
    #                                    subplot_row_col_n=[2,1,1],
    #                                    gridlines=[np.arange(-51,-10,2),
    #                                               np.arange(100,150,4)])
    # add coastline
    ax.coastlines()
    plt.title("Waroona synoptic")
    
    ## Add box around zoomed in area
    xy = [inner[0], inner[2]]
    width = inner[1]-inner[0]
    height = inner[3]-inner[2]
    ax.add_patch(patches.Rectangle(xy=xy,
                                   width=width,
                                   height=height,
                                   #facecolor=None,
                                   fill=False,
                                   edgecolor='blue',
                                   linewidth=2,
                                   #linestyle='-',
                                   alpha=0.6, 
                                   transform=proj))
    ## add text?
    
    ## Add scale
    scaleloc=(0.2,0.05)
    plotting.scale_bar(ax,proj,100, location=scaleloc)
    
    ## Look at waroona and yarloop
    #_,ax2,gproj = plotting.map_google(inner, zoom=10, fig=fig,
    #                                  subplot_row_col_n=[2,2,3], draw_gridlines=False)
    _, ax2, proj2 = plotting.map_tiff(locname='waroona',
                                      extent=inner, fig=fig,
                                      subplot_row_col_n=[2,2,3],
                                      show_grid=False, add_locations=True)
    plt.title("Fire location")
    
    ## Add scale
    plotting.scale_bar(ax2,proj2,10, location=scaleloc)
    
    ## Add contour plot showing topography
    plt.subplot(2,2,4)
    plotting.map_topography(inner, topog, latt, lont)
    
    fio.save_fig_to_path("figures/waroona_outline.png",plt,dpi=300)


def show_nests(model_run='waroona_run1', annotate_res=True, title=''):
    """
    show nested grids on stock image of australia, resolution annotated
    """
    #f, ax, gproj = plotting.map_google(plotting._extents_['waroona'], zoom=11)
    locname=model_run.split('_')[0]
    nest=__NESTS__[model_run]
    
    # look at nests within wider view:
    aust = nest['wider_view']
    
    # create figure and projection
    fig, ax, proj = plotting.map_tiff(locname=locname,
                                      extent=aust,
                                      show_grid=True, 
                                      add_locations=False)
    
    # Where are we looking
    #ax.set_extent(aust)
    # add coastline, stock img (tiff is currently just for zoomed in stuff)
    ax.coastlines()
    ax.stock_img()
    
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
    
    outline_waroona()
    
    for mr in ['sirivan_run1','waroona_run1']:
        show_nests(mr)