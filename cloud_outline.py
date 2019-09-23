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
#import matplotlib.colors as col
import matplotlib.pyplot as plt
#import matplotlib.ticker as tick
#import matplotlib.patches as mpatches
import numpy as np
from datetime import datetime,timedelta
# ignore some warnings
import warnings

# local modules
from utilities import plotting, utils, fio, constants

###
## GLOBALS
###
_sn_ = 'cloud_outline'


def clouds_2panel(topog,s,u,v,
                  qc,w,
                  z,lat,lon,
                  dtime,
                  ff = None,
                  extentname='waroona',
                  transect=1, 
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
        cubes need to be passed in
        extentname = { 'waroona' | 'sirivan' } choice of two fire locations
        transect = int from 1 to 6 for transect choice (see figures/transects.png)
        vectorskip reduces quiver density (may need to fiddle)
        quiverscale changes how long the arrows are (also may need to fiddle)
        dtime is datetime of output 
        ext is the plot extension { '.png' | '.eps' }
    '''
    
    ## Plot data, inputs will be [[z],lat,lon]
    
    # datetime timestamp for title
    stitle = dtime.strftime("%Y %b %d %H:%M (UTC)")
    
    # set font sizes
    plotting.init_plots()
    
    # get plot extent, and transect
    extent = plotting._extents_[extentname]
    start,end = plotting._transects_["%s%d"%(extentname,transect)]
    
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
        plt.contour(lon,lat,np.transpose(ff),np.array([0]), 
                    colors='red')
    
    ## Second row is transect plots
    plt.subplot(3,1,2)
    #plotting.transect(theta,z,lat,lon,start,end,topog=topog,
    #                  cmap='plasma',
    #                  contours=np.linspace(290,400,111),
    #                  ztop=ztop)
    wslice, xslice, zslice = plotting.transect_w(w,z,lat,lon,start,end,topog=topog,
                                                 npoints=100,lines=None,
                                                 ztop=ztop)
    plt.title("Vertical motion (m/s)")
    
    qcslice = utils.cross_section(qc,lat,lon,start,end, npoints=100)
    with warnings.catch_warnings():
        # ignore warning when there are no clouds:
        warnings.simplefilter('ignore')
        plt.contour(xslice,zslice,qcslice,np.array([cloud_threshold]),colors='teal')
    
    plt.ylabel('height (m)')
    plt.xlabel('')
    
    plt.subplot(3,1,3)
    plotting.transect_qc(qc,z,lat,lon,start,end,topog=topog,
                        ztop=ztop,)
    # Show transect start and end
    xticks,xlabels = plotting.transect_ticks_labels(start,end)
    plt.xticks(xticks[0::2],xlabels[0::2]) # show transect start and end
    
    # Save figure into animation folder with numeric identifier
    plt.suptitle(stitle)
    
    return plt


def waroona_cloud_loop(dtime,model_version="waroona_run1"):
    '''
    make an hours worth of clouds_2panel plots starting at argument dtime
    First get iris cubes from each of the data files we will read,
        subset the data before reading it to save ram and run faster
        also read fire front at matching times
    then send all the data to plotting method
    '''
    extentname=model_version.split("_")[0]
    extent=plotting._extents_[extentname]
    
    # Read the cubes
    cubes = fio.read_model_run(model_version, fdtime=dtime, 
                               extent=extent,
                               add_z=True,
                               add_winds=True,
                               add_theta=True)
    print("DEBUG:",cubes)
    zth, topog, u,v,s = cubes.extract(['z_th','surface_altitude','u','v','s'])
    qc, w, sh = cubes.extract(['qc','upward_air_velocity','specific_humidity'])
    theta, Ta = cubes.extract(['air_temperature','potential_temperature'])
    cubetimes = utils.dates_from_iris(u)
    lat,lon = w.coord('latitude').points, w.coord('longitude').points
    
    ## fire front
    # read 6 time steps:
    ff = None
    if model_version in ['waroona_run1']:
        ff1, = fio.read_fire(dtimes=cubetimes, extent=extent, firefront=True)
    if model_version in ['waroona_old']:
        zth = zth[0] # old run has time dim as z_th is estimated from p
    
    # loop over different transects
    for i_transect in np.arange(1,6.5,1, dtype=int):
        for tstep in range(len(cubetimes)):
            if model_version in ['waroona_run1']:
                ff = ff1[tstep].data.data
            plt=clouds_2panel(topog.data.data, s[tstep,0].data.data, u[tstep,0].data.data, v[tstep,0].data.data,
                              qc[tstep].data.data, w[tstep].data.data,
                              zth.data.data, lat, lon,
                              dtime=cubetimes[tstep],
                              ff = ff,
                              extentname=extentname,
                              transect=i_transect)
            # figure name and location
            #pname="figures/%s/clouds_outline_X%d/fig_%s%s"%(extentname,transect,dstamp,ext)
            #print("INFO: Saving figure:",pname)
            #plt.savefig(pname,dpi=dpi)
            #plt.close()
            fio.save_fig(model_version,_sn_,cubetimes[tstep],plt,
                         subdir='transect_%d'%i_transect)
            
if __name__ == '__main__':
    
    print("INFO: testing cloud_outline.py")
    #waroona_cloud_loop(datetime(2016,1,5,15))
    for dtime in [ datetime(2016,1,6,6) + timedelta(hours=x) for x in range(2) ]:
        for mv in ['waroona_old','waroona_run1']:
            waroona_cloud_loop(dtime,model_version=mv)

