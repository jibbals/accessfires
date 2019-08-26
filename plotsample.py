#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 13:49:40 2019

  Script to show some maps and cross sections of Waroona
  Will eventually update to be generic and callable with different extents...

  HISTORY:
    8/8/19: created
      jwg574 just read and plot some simple stuff from converted pc file

@author: jesse
"""


# Copy cross section plotting from jeff's code

# FIRST I WANT A MAP 
import matplotlib as mpl
import matplotlib.colors as col
import matplotlib.pyplot as plt
import numpy as np
import iris

# mapping libraries
import cartopy
import cartopy.feature as cpf
# geometry to add onto cartopy
import shapely.geometry as sgeom
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from matplotlib.transforms import offset_copy

from utilities import fio, utils, plotting



def sensible_heat_flux(shfpath = "C:\\Users\\jgreensl\\Desktop\\data\\waroona_fire\\sensible_heat.CSIRO_24h.20160105T1500Z.nc"):
    '''
        check sensible heat flux fire output
    '''
    pname="figures/waroona/heat_flux_1.png"
    print("INFO: reading fire file ",shfpath, '...')
    # read the cube
    shf, = fio.read_nc_iris(shfpath)
    constr_lons = iris.Constraint(longitude = lambda cell: West <= cell <= East )
    constr_lats = iris.Constraint(latitude = lambda cell: South <= cell <= North )
    shf_waroona=shf.extract(constr_lats & constr_lons)
    
    for i in range(4):
        plt.subplot(221+i)
        # add 0.1 so we don't take log of zero (log scale)
        qplt.contourf(shf_waroona[i*450]+0.1, locator=tick.LogLocator())
        plt.gca().coastlines()
        plt.title(['H0', 'H7.5','H15','H22.5'][i])
    plt.suptitle("Heat Flux on day 1")
    plt.savefig(pname)
    print("INFO: Saved figure: ",pname)


