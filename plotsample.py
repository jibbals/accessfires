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
from netCDF4 import Dataset
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


# Plot defaults
mpl.rcParams['font.size'] = 18.0
mpl.rcParams["text.usetex"]      = False     #
mpl.rcParams["legend.numpoints"] = 1         # one point for marker legends
mpl.rcParams["figure.figsize"]   = (12, 10)  #
mpl.rcParams["font.size"]        = 18        # font sizes:
mpl.rcParams["axes.titlesize"]   = 26        # title font size
mpl.rcParams["axes.labelsize"]   = 20        #
mpl.rcParams["xtick.labelsize"]  = 16        #
mpl.rcParams["ytick.labelsize"]  = 16        #
mpl.rcParams['image.cmap'] = 'plasma' #'PuRd' #'inferno_r'       # Colormap default
mpl.rcParams['axes.formatter.useoffset'] = False

pafile = 'data/umnsaa_pa2016010515.nc'
pcfile = 'data/umnsaa_pc2016010515.nc'

# limit for z axis plotting
zl = 4000
# Extents: EWSN
extents = {}
extents['waroona'] = [115.5,116.3, -33.1,-32.5]
extents['waroonas'] = [109,125,-37,-29] # synoptic
latlons = {}
latlons['waroona'] = -32.8430, 115.8526 # latlon of waroona


# ends of cross section
lona = 115.7
lata = -32.8
lonb = 116.3
latb = -32.8
lonc = 115.95
latc = -32.6
lond = 115.85
latd = -33.2

dlat = 0.05
nsec = 5


## Read the model data
#with Dataset(datadir3+'umnsaa_pc201601{:s}.nc'.format(ddhh2),'r') as ncfile: # Jesse's local copy from nci datadir2

with Dataset(pcfile,'r') as ncfile:
    
    ## PULL OUT MASKED ARRAYS
    #zth  = ncfile.variables['height_theta'][0,:,::-1,:] # flip latitudes for interpolate
    #zrho = ncfile.variables['height_rho'  ][0,:,::-1,:]
    lat  = ncfile.variables['latitude'   ][:]
    lon  = ncfile.variables['longitude'  ][:]
    lat1 = ncfile.variables['latitude_0' ][:] # Edges??
    lon1 = ncfile.variables['longitude_0'][:] # Edges??
    z   = ncfile.variables['model_level_number'  ][:]
    z1  = ncfile.variables['model_level_number_0'][:]
    Ta  = ncfile.variables['air_temperature'][:,:] # at surface? only stored at surface?
    p   = ncfile.variables['air_pressure'][0,0:70,:,:] # in pascals
    pmsl = ncfile.variables['air_pressure_at_sea_level'][0,:,:]
    u1  = ncfile.variables['x_wind'][0,0:70,:,:] # wind speeds are on their directional grid edges (z,lats,lon1)
    v1  = ncfile.variables['y_wind'][0,0:70,:,:] # [z,lat1,lons]
    q   = ncfile.variables['specific_humidity_0'][0,0:70,:,:]
    #w   = ncfile.variables['upward_air_velocity'][0,0:70,:,:]
    #qc  = ncfile.variables['mass_fraction_of_cloud_liquid_water_in_air'][0,0:70,:,:] + ncfile.variables['mass_fraction_of_cloud_ice_in_air'][0,0:70,:,:]

#nz,ny,nx = Ta.shape
nz,ny,nx = p.shape

# READ TOPOG DATA FROM PA
with Dataset(pafile,'r') as ncfile:
    topog = ncfile.variables['surface_altitude'][:,:]
    latt = ncfile.variables['latitude' ][:]
    lont = ncfile.variables['longitude'][:]


## Destagger winds
#u = np.tile(np.nan,(nz,ny,nx))
u = np.tile(np.nan,(nz,ny,nx)) # tile repeats the nan accross nz,ny,nx dimensions
u[:,:,1:] = 0.5*(u1[:,:,1:] + u1[:,:,:-1]) # interpolation of edges
v = 0.5*(v1[:,1::,] + v1[:,:-1,:]) # interpolation of edges
s = np.hypot(u,v) # Speed is hypotenuse of u and v
lonu = lon1
latu = lat
lonv = lon
latv = lat1

# Fudge some height data
zth = -(287*300/9.8)*np.log(p/pmsl[np.newaxis,:,:])
zrho = zth

theta = Ta*(1e5/p)**(287.05/1004.64)

thcon = np.arange(280,320,2) # Temperature axis?
cmapth = plt.cm.get_cmap('YlOrRd')

wcon = 2.0**np.arange(-2,6)
wcon = np.union1d(np.union1d(wcon,-wcon),np.array([0]))
wcon2 = np.array([0])
cmapw = plt.cm.get_cmap('PiYG')
cmapw.set_over('k')
wnorm = col.SymLogNorm(0.25)  # symmetrical logarithmic scale is logarithmic in both the positive and negative directions from the origin

scon = np.arange(0,22,2)
cmaps = plt.cm.get_cmap('YlGnBu')

qccon = np.arange(0.0,2.25,0.25)
cmapqc = plt.cm.get_cmap('YlGnBu')

## First plot showing area with grid and site

proj = cartopy.crs.PlateCarree()
zoomedland = cartopy.feature.NaturalEarthFeature(category='physical',
                                                 name='land',
                                                 scale='10m',
                                                 facecolor="none")
zoomedcoast = cartopy.feature.NaturalEarthFeature(category='physical',
                                                  name='coastline',
                                                  scale='10m')
states= cartopy.feature.NaturalEarthFeature(
          category='cultural',
          name='admin_1_states_provinces_lines',
          scale='50m',
          facecolor='none')


plt.close()
# Create a Stamen terrain background instance.
stamen_terrain = cimgt.Stamen('terrain-background')
fig = plt.figure()
# Create a GeoAxes in the tile's projection.
ax = fig.add_subplot(1, 1, 1, projection=proj)#stamen_terrain.crs)

# Limit the extent of the map to a small longitude/latitude range.
ax.set_extent(extents['waroona'], crs=proj)

# Add the Stamen data at zoom level 14.
ax.add_image(stamen_terrain, 11)

# Grid lines and labels
gl = ax.gridlines(crs=proj, linewidth=1, color='black', alpha=0.5, linestyle='--', draw_labels=True)
gl.xlabels_top = False
gl.ylabels_left = False
gl.ylabels_right=True
gl.xlines = True
gl.xlocator = mpl.ticker.FixedLocator(np.arange(110,120,0.25))
gl.ylocator = mpl.ticker.FixedLocator(np.arange(-40,-10,0.25))
#gl.xformatter = LONGITUDE_FORMATTER
#gl.yformatter = LATITUDE_FORMATTER
#gl.xlabel_style = {'color': 'red', 'weight': 'bold'}


# Add a marker for waroona volcano. LON, LAT, ... argument order important!
plt.plot(latlons['waroona'][1], latlons['waroona'][0],
         color='red', linewidth=0, marker='o',
         transform=proj)
dx,dy = 0.025,0.025
plt.text(latlons['waroona'][1]+dx, latlons['waroona'][0]+dy, 'Waroona',
         horizontalalignment='right',
         transform=proj)

# Show synoptic area for reference
sub_ax = plt.axes([0.65, 0.65, 0.2, 0.2], projection=proj)
sub_ax.set_extent(extents['waroonas'], crs=proj)
# Add the land, coastlines.
#sub_ax.add_feature(zoomedland)
# Add the Stamen data.
sub_ax.add_image(stamen_terrain, 7)
sub_ax.add_feature(states,zorder=2)
#sub_ax.add_feature(zoomedcoast)
# sgeom boxes based on ESWN?
drawbox=extents['waroona']
extent_box = sgeom.box(drawbox[0], drawbox[2], drawbox[1], drawbox[3])
sub_ax.add_geometries([extent_box], proj, color='none',
                          edgecolor='blue', linewidth=2)


## finally show or save figure
##

plt.show()
