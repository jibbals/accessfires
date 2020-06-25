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
from datetime import datetime
from scipy import interpolate

# local modules
from utilities import plotting, utils, fio

###
## GLOBALS
###
_sn_ = 'emberstorm'

# transects for emberstorm plots: [lat0,lon0, lat1,lon1]
x0,x1 = 115.87, 116.02
y0=-32.87
yspread = np.array([0, .008, -.008, .02, -.02, .032,-.032])
_transects_ = [ [[y,x0],[y,x1]] for y in y0+yspread ]

###
## METHODS
###


def emberstorm(theta, u, v, w, z, topog,
               lat, lon,
               ff=None,
               wmap=None,
               wmap_height=None,
               ztop=700,
               transect=_transects_[0],
               shadows=None):
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
        wmap: array[lat,lon] - optional additional vertical motion map to be put on top of topography
        transect: [lat0,lon0],[lat1,lon1] - for transect that will be shown
        shadows: [transects] - faint transects to add to topography (for style)
    '''
    
    extent = [lon[0],lon[-1],lat[0],lat[-1]]
    
    ## figure
    fig, axes = plt.subplots(3,1)
    plt.sca(axes[0])
    ## First plot, topography
    cs,cb = plotting.map_topography(extent,topog,lat,lon,title="Overview")
    
    start, end = transect
    
    # Add fire front contour
    ff_lead_frac=None
    if ff is not None:
        plotting.map_fire(ff,lat,lon, transform=False)
        # west most burnt bit is fire front lead
        # grab first index from left to right along burnt longitudes
        if np.sum(ff<=0) > 0:
            #print("debug:", ff.shape, len(lat), len(lon), np.sum(ff<=0,axis=0).shape)
            ff_lead = np.where(np.sum(ff<=0, axis=1)>0)[0][0]
            print("debug: fire front most western point is ", lon[ff_lead])
            ff_lead_frac=(start[1]-lon[ff_lead])/(start[1] - end[1])
            print("debug: %.2f%% of the way along the transect"%ff_lead_frac)
            
    
    # start to end x=[lon0,lon1], y=[lat0, lat1]
    plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
             linewidth=2,) 
             #marker='>', markersize=7, markerfacecolor='white')
    
    # Add markers where other transects will be 
    if shadows is not None:
        for sstart,send in shadows:
            plt.plot([sstart[1],send[1]],[sstart[0],send[0], ], 
                     color='b', alpha=0.4, linestyle='None',
                     marker='>', markersize=2)
    
    # add nearby towns
    plotting.map_add_locations(['waroona'], text=[''], color='grey',
                               marker='o', markersize=4)
    plotting.map_add_locations(['fire_waroona'], text=[''], 
                               marker='*', color='r')
    
    # show winds
    plt.streamplot(lon,lat,u[0],v[0], color='k',
                   density=(.8,.5),
                   )#alpha=0.7)
    
    
    # Finally add some faint pink or green hatching based on vertical wind motion
    if wmap is not None:
        with warnings.catch_warnings():
            # ignore warning when there are no contours to plot:
            warnings.simplefilter('ignore')
            plt.contour(lon, lat, -1*wmap, levels=[1,3], 
                        linestyles=['dashed','solid'],
                        colors=('aquamarine',))
            plt.contour(lon, lat, wmap, levels=[1,3], 
                        linestyles=['dashed','solid'],
                        colors=('pink',))
            #plt.contour(lon,lat,wmap, levels=[-2.1,2.1],
            #            alpha=0.5,
            #            colors=['aquamarine','pink'])
            axes[0].annotate('vertical motion at %dm altitude'%wmap_height, 
                xy=[0.01,-.05], xycoords='axes fraction', fontsize=8)
            axes[0].annotate('dashed is 1m/s, solid is 3m/s', 
                xy=[0.01,-.11], xycoords='axes fraction', fontsize=8)
    
    # cut down to desired extent
    plt.ylim(extent[2:]); plt.xlim(extent[0:2])
    
    ## Turn off the tick values
    plt.xticks([]); plt.yticks([])
    ## Subplot 2, transect of potential temp
    # only looking up to 1km
    # horizontal points
    npoints = utils.number_of_interp_points(lat,lon,start,end)
    plt.sca(axes[1])
    trets = plotting.transect_theta(theta, z, lat, lon, start, end, npoints=npoints,
                                    topog=topog, ff=ff, ztop=ztop,
                                    contours=np.arange(290,320.1,0.5),
                                    lines=None, #np.arange(290,321,2), 
                                    linestyles='dashed')
    # add faint lines for clarity
    thetaslice,xslice,zslice=trets
    plt.contour(xslice,zslice,thetaslice,np.arange(290,320.1,1),colors='k',
                alpha=0.5, linestyles='dashed', linewidths=0.5)
    if wmap_height is not None:
        axes[1].annotate('', xy=(xslice[0,0],wmap_height), 
            xytext=(xslice[0,5], wmap_height),
            arrowprops=dict(facecolor='grey', arrowstyle='wedge,tail_width=0.5', alpha=0.5),
            fontsize=8)
    
    ## Next plot
    plt.sca(axes[2])
    _,slicex,slicez = plotting.transect_w(w,z,lat,lon,start,end,
                                          npoints=npoints,topog=topog, ff=ff,
                                          ztop=ztop)
    # add east-west-vertical winds
    uslice,_,_ = utils.transect(u,lat,lon,start,end,nx=npoints)
    #vslice = utils.cross_section(v,lat,lon,start,end,npoints=npoints)
    wslice,_,_ = utils.transect(w,lat,lon,start,end,nx=npoints)
    
    ## Streamplot of east-west and vertical winds
    plotting.streamplot_regridded(slicex,slicez,uslice,wslice,
                                  density=(.8,.5))
    
    # reset plot edges?
    #plt.ylim(extent[2:]); plt.xlim(extent[0:2])
    
    if ff_lead_frac is not None:
        if ff_lead_frac > 1: ff_lead_frac=1
        if ff_lead_frac < 0: ff_lead_frac=0
        for ax in [axes[1],axes[2]]:
            ax.annotate('', xy=[ff_lead_frac,0.0],
                        xytext=(ff_lead_frac,0.09),
                        arrowprops=dict(facecolor='red', arrowstyle='wedge,tail_width=0.5', alpha=0.5),
                        xycoords='axes fraction', fontsize=8, color='red',
                        zorder=10)
    
    return fig, axes

def topdown_emberstorm(fig=None, subplot_row_col_n=None, 
                       extent=[115.71, 116.1, -33.05,-32.78],
                       lats=None,lons=None, ff=None, sh=None, u10=None, v10=None):
    """
    Top down view of Waroona/Yarloop, adding fire front and heat flux and 10m winds
    ARGUMENTS:
        lats/lons are 1D arrays, required if adding other stuff
        ff [lats,lons] : firefront array
        sh [lats,lons] : sensible heat flux
        u10 [lats,lons] : 10m altitude longitudinal winds
        v10 [lats,lons] : 10m altitude latitudinal winds
    RETURNS:
        fig, ax, proj : proj is map projection used by tiff
    """
    # first create map from tiff file
    
    fig, ax, proj = plotting.map_tiff_qgis(
        fname="waroonaz_osm.tiff", 
        extent=extent,
        fig=fig,
        subplot_row_col_n=subplot_row_col_n,
        #EPSG=4326,
        )
    
    if ff is not None:
        # add firefront
        cs_ff = plotting.map_fire(ff,lats,lons,transform=True)
    if sh is not None:
        # add hot spots for heat flux
        cs_sh, cb_sh = plotting.map_sensibleheat(sh,lats,lons)
    
    return fig, ax, proj

def make_plots_emberstorm(model_run='waroona_run1', hours=None, wmap_height=300):
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
            for transecti, transect in enumerate(_transects_):
                # other transects that will be plotted
                shadows = [tran for trani, tran in enumerate(_transects_) if trani!=transecti]
                # current time
                utcstamp = dtime.strftime("%b %d %H:%M (UTC)")
                # winds
                u,v,w = uvw[0][i], uvw[1][i], uvw[2][i]
                # fire
                ffi=None
                if ff is not None:
                    ffi = ff[i].data 
                
                # extra vert map at ~ 300m altitude
                levh  = w.coord('level_height').points
                levhind = np.sum(levh<wmap_height)
                emberstorm(theta[i].data, u.data, v.data, w.data, z.data,
                           topog.data, lat, lon, ff=ffi, 
                           wmap=w[levhind].data, wmap_height=wmap_height,
                           transect=transect, shadows=shadows)
                
                # Save figure into folder with numeric identifier
                stitle="Emberstorm %s"%utcstamp
                plt.suptitle(stitle)
                distance=utils.distance_between_points(transect[0],transect[1])
                plt.xlabel("%.3fE - %.3fE = %.1fkm (at %.3fS)"%(transect[0][1], 
                           transect[1][1], distance/1e3, -1*transect[0][0]))
                fio.save_fig(model_run=model_run, script_name=_sn_, pname=dtime, 
                             plt=plt, subdir='transect%d'%transecti)

if __name__ == '__main__':
    plotting.init_plots()
    mr = 'waroona_run1'
    
    
    hours=fio.model_outputs[mr]['filedates']
    testhours = [datetime(2016,1,6,13)]
    interesting_hours=[datetime(2016,1,6,x) for x in range(8,15)]

    if False:
        # This makes the combined 3 row plot with top down winds and 
        # transects of theta and wind
        make_plots_emberstorm(mr, hours=testhours)
    
    if True:
        # newer plots showing 1: fire + town + winds (based on top panel in make_plots_emberstorm)
        extent=[115.71, 116.1, -33.05,-32.78]
        # read fire
        ff,sh = fio.read_fire(model_run=mr,
                              dtimes=testhours, 
                              extent=extent,
                              sensibleheat=True)
        FF = ff[0].data.data
        SH = sh[0].data.data
        lats = ff.coord('latitude').points
        lons = ff.coord('longitude').points
        
        fig,ax,proj = topdown_emberstorm(extent=extent,
                                         lats=lats,lons=lons,
                                         ff=FF, sh=SH)
        fio.save_fig('test',_sn_,'emberstorm_topdown.png',plt)
