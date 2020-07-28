#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""


# IMPORTS

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches, colors, cm
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection
from datetime import datetime, timedelta
import iris # so we can add level number constraint
import iris.analysis
import iris.quickplot as qplt

# read tiff file
from osgeo import gdal, osr
import cartopy.crs as ccrs

from utilities import utils, plotting, fio, constants

from scipy import interpolate

fdtime = datetime(2016,1,6,12)
fdtime2 = datetime(2016,1,7,4)
mr='waroona_run3'

#cubes=fio.read_model_run(mr,[fdtime,fdtime2], add_z=True,add_winds=True)
#print (cubes)
#z, w, u, topog = cubes.extract(['z_th','upward_air_velocity','u','surface_altitude'])
#lats,lons = w.coord('latitude').points, w.coord('longitude').points
#print(z.shape)

# Read firefront
extent = plotting._extents_['waroonaf'] # synoptic extent
ff, = fio.read_fire(model_run=mr, dtimes=None,
                    extent=extent, firefront=True,
                    HSkip=5)


from shapely import geometry
#import matplotlib.pyplot as plt
#from matplotlib import cm
#import numpy as np
import fiona
import os,json
#from descartes.patch import PolygonPatch

# create some test data with multiple peaks
lon = np.linspace(0,45,100)
lat = np.linspace(-20,32,90)
long,latg=np.meshgrid(lon,lat)
C1=np.sqrt((long-5.)**2+(latg-25)**2)/30.
C2=np.sqrt((long-30.)**2+(latg-1)**2)/10.
m = 30*np.exp(-C1**2)+20.*np.exp(-C2**2)

# make the contourf plot, storing the resulting ContourSet in cs
plt.figure(figsize=[10,5])
plt.subplot(1,2,1)
Nlevels=10
cs = plt.contourf(lon,lat,m,Nlevels,cmap='gist_heat')
plt.title('contourf figure with Nlevels='+str(Nlevels))

# create lookup table for levels
lvl_lookup = dict(zip(cs.collections, cs.levels))

# loop over collections (and polygons in each collection), store in list for fiona
PolyList=[]
for col in cs.collections:
    z=lvl_lookup[col] # the value of this level
    for contour_path in col.get_paths():
        # create the polygon for this level
        for ncp,cp in enumerate(contour_path.to_polygons()):
            lons = cp[:,0]
            lats = cp[:,1]
            new_shape = geometry.Polygon([(i[0], i[1]) for i in zip(lons,lats)])            
            if ncp == 0:
                poly = new_shape # first shape
            else:
                poly = poly.difference(new_shape) # Remove the holes
        PolyList.append({'poly':poly,'props':{'z': z}})

## write the fiona collection

# clean up directories
outname=os.path.join('..','data','shaped_contour')
if os.path.isdir(outname):
    for file in os.listdir(outname):
        os.remove(os.path.join(outname,file))
    os.rmdir(outname)
os.mkdir(outname)

# define ESRI schema, write each polygon to the file
outfi=os.path.join(outname,'shaped_contour.shp')
schema = {'geometry': 'Polygon','properties': {'z': 'float'}}
with fiona.collection(outfi, "w", "ESRI Shapefile", schema) as output:
    for p in PolyList:
        output.write({'properties': p['props'],
            'geometry': geometry.mapping(p['poly'])})

# save the levels and global min/max as a separate json for convenience
Lvls={'levels':cs.levels.tolist(),'min':m.min(),'max':m.max()}
with open(os.path.join(outname,'levels.json'), 'w') as fp:
    json.dump(Lvls, fp)

## Plotting the results: reads data back in, plots the polygons with data only
## from shapefile and levels.txt

ax=plt.subplot(1,2,2)

# read in levels, define colormap
with open(os.path.join(outname,'levels.json')) as jfile:
    Lvls=json.load(jfile)
levels=np.array(Lvls['levels'])
cmap=plt.cm.gist_heat
lv_range=[Lvls['min'],Lvls['max']]

## loop over each shape, pull out level value ('z'), plot a polygon with a color
## matching the colormap.
#ishp=0
#with fiona.open(outfi) as shape:
#    for shp in shape:
#
#        # pull this shape's level and set the color from the map
#        lv=shp['properties']['z'] # this shape's level
#        clr=cmap((lv - lv_range[0])/(lv_range[1]-lv_range[0]))
#
#        # build the polygon and add the patch
#        coords=shp['geometry']['coordinates'][0] # coords of this polygon
#        poly=geometry.Polygon(coords)
#        patch = PolygonPatch(poly, facecolor=clr, edgecolor=clr)
#        ax.add_patch(patch)
#
#        # track max/min coordinate values
#        bnds=poly.bounds
#        rng_C={'lon':{'min':bnds[0],'max':bnds[2]},
#               'lat':{'min':bnds[1],'max':bnds[3]}}
#        if ishp==0:
#            rngs=rng_C
#        else:
#            for ll in ['lon','lat']:
#                rngs[ll]['max']=max([rngs[ll]['max'],rng_C[ll]['max']])
#                rngs[ll]['min']=min([rngs[ll]['min'],rng_C[ll]['min']])
#        ishp=ishp+1
#
#ax.set_xlim([rngs['lon']['min'],rngs['lon']['max']])
#ax.set_ylim([rngs['lat']['min'],rngs['lat']['max']])
#plt.title('polygon patches from shapefile')

#plt.show()