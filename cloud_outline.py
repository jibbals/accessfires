# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 15:58:21 2019
    Create cloud summary plots
@author: jgreensl
"""

import matplotlib
# don't plot on screen, send straight to file
# this is for lack of NCI display
matplotlib.use('Agg', warn=False)

# plotting stuff
from matplotlib import colors
import matplotlib.pyplot as plt
#import matplotlib.ticker as tick
#import matplotlib.patches as mpatches
import numpy as np
from datetime import datetime

# ignore some warnings
import warnings

# local modules
from utilities import plotting, utils, fio, constants

###
## GLOBALS
###
_sn_ = 'cloud_outline'


def clouds_2panel(topog,s,u,v,
                  qc,th,
                  z,lat,lon,
                  transect,
                  ff = None,
                  extentname='waroona',
                  nquivers=20,
                  quiverscale=60,
                  ztop=10000,
                  ext='.png',
                  cloud_threshold=constants.cloud_threshold,
                  dpi=200,
                  ):
    '''
    311 Plot showing windspeed map and wind speed, along with near sites and transect
    312 relative humidity along transect
    313 water and ice content along transect
    INPUTS:
        topog: [lat,lon] surface altitude (m)
        s,u,v: [z, lat, lon] wind speed, x-wind, y-wind (m/s)
        qc,th: [z,lat,lon] cloud content, potential temp (g/kg, K)
        z,lat,lon: dimensions for other inputs
        extentname = { 'waroona' | 'sirivan' } choice of two fire locations
        transect = int from 1 to 6 for transect choice (see figures/transects.png)
        vectorskip reduces quiver density (may need to fiddle)
        quiverscale changes how long the arrows are (also may need to fiddle)
        dtime is datetime of output 
        ext is the plot extension { '.png' | '.eps' }
    '''
    
    # set font sizes
    plotting.init_plots()
    
    # get plot extent, and transect
    extent = [lon[0],lon[-1],lat[0],lat[-1]]
    #extent = plotting._extents_[extentname]
    start,end = transect
    
    plt.figure(figsize=[7,10])
    plt.subplot(3,1,1)
    
    # top panel is wind speed surface values
    plotting.map_contourf(extent,s,lat,lon,
                          clevs = np.linspace(0,15,31),
                          cmap=plotting._cmaps_['windspeed'],
                          clabel='m/s')
    plt.title('Horizontal wind speed')
    
    # start to end x=[lon0,lon1], y=[lat0, lat1]
    plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
             linewidth=2, 
             marker='X', markersize=7,markerfacecolor='white')
    
    # add nearby towns
    plotting.map_add_locations_extent(extentname)
    
    # Add vectors for winds
    # just surface, and one every N points to reduce density
    vsu = len(lon)//nquivers
    vsv = len(lat)//nquivers
    skip = (slice(None,None,vsv),slice(None,None,vsu))
    #mlon,mlat = np.meshgrid(lon,lat)
    plt.quiver(lon[skip[1]],lat[skip[0]],u[skip],v[skip], scale=quiverscale)
    
    # Add fire outline
    if ff is not None:
        with warnings.catch_warnings():
            # ignore warning when there is no fire:
            warnings.simplefilter('ignore')
            plt.contour(lon,lat,np.transpose(ff),np.array([0]), 
                        colors='red')
    
    ## Second row is transect plots
    plt.subplot(3,1,2)
    #Lets do theta using log normal after 300 degrees
    #norm = colors.SymLogNorm(300)
    #contours = np.arange(280,350,1)
    #lines = np.union1d(np.arange(280,301,2), np.arange(310,351,10))
    thslice, xslice, zslice = plotting.transect_theta(th,z,lat,lon,start,end,topog=topog,
                                                      npoints=100,
                                                      lines=None,
                                                      ztop=ztop)
    plt.title("Vertical motion (m/s)")
    
    qcslice = utils.cross_section(qc,lat,lon,start,end, npoints=100)
    with warnings.catch_warnings():
        # ignore warning when there are no clouds:
        warnings.simplefilter('ignore')
        plt.contour(xslice,zslice,qcslice,np.array([cloud_threshold]),colors='teal')
    
    plt.ylabel('height (m)')
    plt.xlabel('')
    
    ## subplot with cloud contourf
    plt.subplot(3,1,3)
    plotting.transect_qc(qc,z,lat,lon,start,end,topog=topog,
                        ztop=ztop,)
    # Show transect start and end
    xticks,xlabels = plotting.transect_ticks_labels(start,end)
    plt.xticks(xticks[0::2],xlabels[0::2]) 


def cloud_outline_model(model_run = 'waroona_run1', dtime=datetime(2016,1,5,15)):
    '''
    make an hours worth of clouds_2panel plots starting at argument dtime
    Read model run, send good bits to be plotted over our 6 transects
    '''
    extentname=model_run.split("_")[0]
    extent=plotting._extents_[extentname]
    
    # Read the cubes
    cubes = fio.read_model_run(model_run, fdtime=dtime, 
                               extent=extent,
                               add_z=True,
                               add_winds=True,
                               add_theta=True)
    
    zth, topog, u,v,s = cubes.extract(['z_th','surface_altitude','u','v','s'])
    qc, w, sh = cubes.extract(['qc','upward_air_velocity','specific_humidity'])
    T, theta = cubes.extract(['air_temperature','potential_temperature'])
    cubetimes = utils.dates_from_iris(u)
    lat,lon = w.coord('latitude').points, w.coord('longitude').points
    
    ## fire front
    # read 6 time steps:
    ff, = fio.read_fire(model_run, dtimes=cubetimes, extent=extent, firefront=True)
    
    # loop over available timesteps
    topogi = topog.data.data
    zi = zth.data.data
    for i,dt in enumerate(cubetimes):
        # datetime timestamp for title
        stitle = dt.strftime("%Y %b %d %H:%M (UTC)")
        si, ui, vi, qci, thi, ffi = [s[i,0].data.data, u[i,0].data.data, 
                                        v[i,0].data.data, qc[i].data.data, 
                                        theta[i].data.data, ff[i].data.data]
        # loop over transects
        for ii in range(6):
            transect = plotting._transects_['%s%d'%(extentname,ii+1)]
            clouds_2panel(topogi,si,ui,vi, qci,thi, zi, lat, lon,
                          transect=transect, ff = ffi, extentname=extentname)
            # Save figure into animation folder with numeric identifier
            plt.suptitle(stitle)
            fio.save_fig(model_run, _sn_, dt, plt,
                         subdir='transect_%d'%(ii+1))

if __name__ == '__main__':
    
    testing=False
    
    for mv in ['sirivan_run1', 'waroona_old', 'waroona_run1']:
        datetimes = fio.model_outputs[mv]['filedates']
        if testing:
            datetimes = datetimes[0:2]
        
        for dt in datetimes:
            cloud_outline_model(mv, dt)

