#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 10:09:53 2019
    Examine emberstorm area winds and weather
@author: jesse
"""

import matplotlib
matplotlib.use('Agg',warn=False)

# plotting stuff
import matplotlib.pyplot as plt
import numpy as np

import warnings

# local modules
from utilities import plotting, utils, fio

###
## GLOBALS
###
_sn_ = 'emberstorm'

###
## METHODS
###


def emberstorm(theta, u, v, w, z, topog,
               lat, lon,
               ff=None,
               nquivers=15,
               quiverscale=60,
               ztop=700,
               transect=plotting._transects_['emberstorm1']):
    '''
    311: contours map, overset with surface wind vectors and transect line
    312: potential temperature (theta) cross section
    313: horizontal winds transect, overlaid with vertical motion contours
    INPUTS:
        theta: array[lev,lat,lon] - potential temp (K)
        u,v,w: array[lev,lat,lon] - x,y,z wind speeds
        z: array[lev,lat,lon] - height of model levels
        topog: array[lat,lon] - surface altitude over lat, lon dim
        lat,lon: 1d arrays
        ff: array[lat,lon] - fire front
    '''
    extent = [lon[0],lon[-1],lat[0],lat[-1]]
    
    ## figure
    fig, axes = plt.subplots(3,1)
    plt.sca(axes[0])
    ## First plot, topography
    cs,cb = plotting.map_topography(extent,topog,lat,lon,title="")
    # set colorbar ticks to scalars..
    #tick.ScalarFormatter()
    
    # Add fire front contour
    if ff is not None:
        plotting.map_fire(ff.T,lat,lon)
    
    start, end = transect
    # start to end x=[lon0,lon1], y=[lat0, lat1]
    plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
             linewidth=2, 
             marker='>', markersize=7, markerfacecolor='white')
    
    # add nearby towns
    plotting.map_add_locations(['waroona'], text=['Waroona'], color='grey',
                               marker='o', markersize=4)
    plotting.map_add_locations(['fire_waroona'], text=[''], 
                               marker='*', color='r')

    # Quiver, reduce arrow density
    plotting.map_quiver(u[0],v[0],lat,lon,nquivers=nquivers,
                        pivot='middle', scale=quiverscale)
    
    ## Subplot 2, transect of potential temp
    # only looking up to 1km
    # horizontal points
    npoints = 100
    plt.sca(axes[1])
    trets = plotting.transect_theta(theta, z, lat, lon, start, end, npoints=npoints,
                                    topog=topog, ff=ff, ztop=ztop,
                                    contours=np.arange(290,320.1,0.5),
                                    lines=None, #np.arange(290,321,2), 
                                    linestyles='dashed')
    # add faint lines for visibility
    thetaslice,xslice,zslice=trets
    plt.contour(xslice,zslice,thetaslice,np.arange(290,320.1,1),colors='k',
                alpha=0.5, linestyles='dashed', linewidths=0.5)
    
    ## Next plot
    plt.sca(axes[2])
    _,xslice,zslice = plotting.transect_w(w,z,lat,lon,start,end,
                                          npoints=npoints,topog=topog, ff=ff,
                                          ztop=ztop)
    # add horizontal winds
    uslice = utils.cross_section(u,lat,lon,start,end,npoints=npoints)
    #vslice = utils.cross_section(v,lat,lon,start,end,npoints=npoints)
    wslice = utils.cross_section(w,lat,lon,start,end,npoints=npoints)
    
    # lets cut away the upper levels before making our quiver
    ztopirows, ztopicols = np.where(zslice < ztop) # index rows and columns
    ztopi = np.max(ztopirows) # highest index with height below ztop
    
    # quiver east-west and vertical winds
    plotting.map_quiver(uslice[:ztopi+4,:], wslice[:ztopi+4,:], 
                        zslice[:ztopi+4,:], xslice[:ztopi+4,:], 
                        nquivers=nquivers, scale=quiverscale*3,
                        alpha=0.5, pivot='middle')
    
    return fig, axes
    
def make_plots_emberstorm(model_run='waroona_run1', hours=None):
    """
    run emberstorm plotting method on model output read by my iris fio scripts
    """
    extent = plotting._extents_['waroona']
    # read topog
    topog=fio.read_topog(model_run,extent=extent)
    lat = topog.coord('latitude').points
    lon = topog.coord('longitude').points
    # Read model run
    umdtimes = fio.model_outputs[model_run]['filedates']
    if hours is not None:
        umdtimes=hours
    # read one model file at a time
    for umdtime in umdtimes:
        cubelist = fio.read_model_run(model_run, 
                                      fdtime=umdtime,
                                      extent=extent,
                                      add_topog=False,
                                      add_z=True, add_theta=True, add_winds=True)
        
        theta, = cubelist.extract('potential_temperature')
        dtimes = utils.dates_from_iris(theta)
        # read fire
        ff, = fio.read_fire(model_run=model_run, dtimes=dtimes, extent=extent)
        
        # pull out bits we want
        uvw = cubelist.extract(['u','v','upward_air_velocity'])
        z, = cubelist.extract(['z_th']) # z has no time dim
        # for each time slice pull out potential temp, winds
        for i,dtime in enumerate(dtimes):
            for transecti, transect in enumerate([plotting._transects_['emberstorm%d'%d] for d in range(1,4)]):
                utcstamp = dtime.strftime("%b %H:%M (UTC)")
                u,v,w = uvw[0][i], uvw[1][i], uvw[2][i]
                ffi=None
                if ff is not None:
                    ffi = ff[i].data.T # fire needs to be transposed...
                #print("DEBUG:",[np.shape(arr) for arr in [theta, u, v, w, z, topog, ff, lat, lon]])
                emberstorm(theta[i].data,u.data,v.data,w.data,z.data,topog.data,
                           lat, lon, ff=ffi, transect=transect)
                
                # Save figure into folder with numeric identifier
                stitle="Emberstorm %s"%utcstamp
                plt.suptitle(stitle)
                distance=utils.distance_between_points(transect[0],transect[1])
                plt.xlabel("%.1fkm transect"%(distance/1e3))
                fio.save_fig(model_run=model_run, script_name=_sn_, pname=dtime, 
                             plt=plt, subdir='transect%d'%transecti)

if __name__ == '__main__':
    plotting.init_plots()
    mr = 'waroona_run2'
    hours=fio.model_outputs[mr]['filedates']#[-6:] # can just do last 6 hours
    make_plots_emberstorm(mr, hours=hours)