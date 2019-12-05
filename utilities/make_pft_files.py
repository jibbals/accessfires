# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 14:56:41 2019

    Call this file from main working directory (accessfires)
    This should create PFT.nc data for each model run but will take 45 million 
    years so make sure you're using short interactive queeueuueu
    
@author: jgreensl
"""

import numpy as np

import PFT_work
import iris
from utilities import fio

# deresolved horizontally (don't need high res)
Hskip=5

for mr in ['waroona_run1','waroona_old','sirivan_run1']:
    mr='waroona_run1'
    hours = fio.model_outputs[mr]['filedates']
    PFT_full=[]
    dtimes=[]
    for hour in hours:
        # Read the cubes for one hour at a time
        cubes = fio.read_model_run(mr, fdtime=hour,
                                   add_z=True, add_RH=True,
                                   add_topog=True, add_winds=True,
                                   add_theta=True)
        
        PFT = PFT_work.PFT_from_cubelist(cubes, latskip=Hskip, lonskip=Hskip)
        PFT_full.append(PFT)
        dtimes.append(cubes[0].dim_coords[0].points)
    
    # extend along time dimension...
    PFT=np.concatenate(PFT_full,axis=0)
    dtimes = np.concatenate(dtimes,axis=0) # array of seconds since 1970-01-01 00:00:00
    lats = cubes[0].coord('latitude').points[::Hskip]
    lons = cubes[0].coord('longitude').points[::Hskip]
    Pa = cubes[0][:,0,::Hskip,::Hskip]
    
    # copy dimensions
    timecoord = iris.coords.DimCoord(dtimes, 'time', units=Pa.coord('time').units)
    dim_coords = [(b,a) for [a,b] in enumerate(Pa.dim_coords)]
    # update time dimension
    dim_coords[0] = (timecoord,0)
    # create cube from PFT array
    PFTcube = iris.cube.Cube(PFT,
                             var_name='PFT',
                             units='Gigawatts',
                             dim_coords_and_dims=dim_coords)
    # Keep some attributes too
    attrkeys=['source', 'um_version', 'institution', 'title', 
              'summary', 'project', 'acknowledgment', 
              'license', 
              'geospatial_lon_units',
              'geospatial_lat_units',
              'publisher_institution', 'publisher_name', 'publisher_type', 
              'publisher_url', 'naming_authority']
    attributes = { k:Pa.attributes[k] for k in attrkeys }
    attributes['geospatial_lon_max'] = np.max(lons)
    attributes['geospatial_lat_max'] = np.max(lats)
    PFTcube.attributes=attributes
    
    # save file
    fname='data/%s/PFT.nc'%mr
    iris.save(PFTcube,fname)

    # test file:
    f = iris.load(fname)
    print("INFO: SAVED ",fname)
    print(f)