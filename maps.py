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
from matplotlib import patheffects

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
                             'tiffname':'WesternAus.tiff',
                             'pointnames':['perth',]
                             },
            'sirivan_run1':{'centre':[-32.001400, 149.798600],
                            'resolution':[.036,.01,.0028],
                             'nlats':[384,576,576],
                             'nlons':[384,576,576],
                             'wider_view':[125,162,-44,-11],
                             },
            }
# old and run2 are the same as run1
__NESTS__['waroona_run2'] = __NESTS__['waroona_run1']
__NESTS__['waroona_run3'] = __NESTS__['waroona_run1']
__NESTS__['waroona_old'] = __NESTS__['waroona_run1']

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
    
    cube_topog = fio.read_topog('waroona_run3',inner)
    latt = cube_topog.coord('latitude').points
    lont = cube_topog.coord('longitude').points
    topog = cube_topog.data.data
    
    # Google map image tiles view of synoptic map
    fig, ax, proj = plotting.map_tiff(locname='waroona',
                                      extent=extent, #fig=fig,
                                      subplot_row_col_n=[2,1,1],
                                      show_grid=False, locnames=['waroona','perth'])
    
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
    _, ax2, proj2 = plotting.map_tiff(locname='waroona',
                                      extent=inner, fig=fig,
                                      subplot_row_col_n=[2,2,3],
                                      show_grid=False, locnames=['waroona','yarloop','fire_waroona'])
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
    tiffname = __NESTS__[model_run]['tiffname']
    locnames = __NESTS__[model_run]['pointnames']
    # look at nests within wider view:
    aust = nest['wider_view']
    
    # create figure and projection
    fig, ax, proj = plotting.map_tiff_qgis(tiffname,
                                           extent=aust,
                                           show_grid=True, 
                                           locnames=locnames)
    
    # add coastline, stock img (tiff is currently just for zoomed in stuff)
    #ax.coastlines()
    #ax.stock_img()
    
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
        # transform for our epsg
        maptransform=cartopy.crs.PlateCarree()._as_mpl_transform(ax)
        
        # outline for text/markers
        outlinecolor='k'
        #marker_effects=[patheffects.Stroke(linewidth=5, foreground=outlinecolor), patheffects.Normal()]
        text_effects = [patheffects.withStroke(linewidth=3, foreground=outlinecolor)]
        text_color="wheat"
        ## add text
        txt = ax.annotate(['Nest 1','Nest 2','N3'][i], xy=botleft, 
                          xycoords=maptransform, 
                          color=text_color,
                          ha='left', va='bottom')
        txt.set_path_effects(text_effects)
        
        if annotate_res:
            nestres = [3.5, 1.0, 0.3][i]
            txt = ax.annotate('Nest %d: %0.1fx%0.1f km squares'%(i+1, nestres,nestres),
                              xy=[.055, .95 - .05*i],
                              xycoords='axes fraction',
                              color=text_color, fontsize=13,
                              ha='left', va='top')
            txt.set_path_effects(text_effects)
            
    if title=='':
        locname = str.capitalize(model_run.split('_')[0])
        if locname == 'Sirivan':
            locname = 'Sir Ivan'
        title = "Fire-Access model run nests at %s"%locname
    plt.title(title)
    ## Add zoom of nest 3?
    fio.save_fig(model_run,_sn_,'nested_grid',plt)
    
if __name__=='__main__':
    
    #outline_waroona()
    
    for mr in ['waroona_run3']:
        show_nests(mr)
