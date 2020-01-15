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
               ztop=1000,
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
        with warnings.catch_warnings():
            # ignore warning when there are no fires:
            warnings.simplefilter('ignore')
            plt.contour(lon,lat,np.transpose(ff),np.array([0]), 
                        colors='red')
    
    start, end = transect
    # start to end x=[lon0,lon1], y=[lat0, lat1]
    plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
             linewidth=2, 
             marker='>', markersize=7, markerfacecolor='white')
    
    # add nearby towns
    plotting.map_add_locations(['waroona'],text=['Waroona'], textcolor='k', 
                               marker='o', color='grey', markersize=4)
    plotting.map_add_locations(['fire_waroona'],text=[''], 
                               marker='*', color='r')

    # Quiver, reduce arrow density
    vsu = len(lon)//nquivers
    vsv = len(lat)//nquivers
    skip = (slice(None,None,vsv),slice(None,None,vsu))
    
    # wind vectors quiver
    plt.quiver(lon[skip[1]],lat[skip[0]],u[0][skip],v[0][skip], 
               scale=quiverscale, pivot='middle')
    
    ## Subplot 2, transect of potential temp
    # only looking up to 1km
    # horizontal points
    npoints = 100
    plt.sca(axes[1])
    plotting.transect_theta(theta,z,lat,lon,start,end,npoints=npoints,topog=topog,
                            ztop=ztop,
                            cmap='gist_rainbow',
                            contours=np.arange(290,330.1,0.5),
                            lines=np.arange(290,331,2))
    
    plt.sca(axes[2])
    
    _,xslice,zslice = plotting.transect_w(w,z,lat,lon,start,end,npoints=npoints,topog=topog,
                                          ztop=ztop)
    # add horizontal winds
    uslice = utils.cross_section(u,lat,lon,start,end,npoints=npoints)
    vslice = utils.cross_section(v,lat,lon,start,end,npoints=npoints)
    
    #vskip = slice(None,None,np.shape(z)[0]//14) # just leave around 14 vertical 
    #print("DEBUG: horizontal slice shapes",np.shape(uslice),np.shape(vslice), np.shape(xslice), np.shape(zslice))
    # want twice as many quivers (vertically) for this very zoomed in slice
    nz,_ = np.shape(uslice)
    skip = (slice(None,None,nz//(nquivers*2)),slice(None,None,npoints//nquivers))
    plt.quiver(xslice[skip], zslice[skip], uslice[skip], vslice[skip], 
               scale=quiverscale*4, alpha = 0.5, pivot='middle')
    
    return fig, axes
    
def make_plots_emberstorm(model_run='waroona_run1'):
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
                    ffi = ff[i].data
                #print("DEBUG:",[np.shape(arr) for arr in [theta, u, v, w, z, topog, ff, lat, lon]])
                emberstorm(theta[i].data,u.data,v.data,w.data,z.data,topog.data,
                           lat,lon,ff=ffi, transect=transect)
                
                # Save figure into folder with numeric identifier
                stitle="Emberstorm %s"%utcstamp
                plt.suptitle(stitle)
                fio.save_fig(model_run=model_run, script_name=_sn_, pname=dtime, 
                             plt=plt, subdir='transect%d'%transecti)

if __name__ == '__main__':
    
    print("INFO: testing cloud_outline.py")
    
    [ make_plots_emberstorm(mr) for mr in ['waroona_run2','waroona_run1','waroona_old'] ]