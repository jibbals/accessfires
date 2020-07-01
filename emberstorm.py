#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 10:09:53 2019
    Examine emberstorm area winds and weather
@author: jesse
"""

import matplotlib
#matplotlib.use('Agg',warn=False)

# plotting stuff
import matplotlib.pyplot as plt
import numpy as np
import warnings
from datetime import datetime,timedelta
from scipy import interpolate
import cartopy.crs as ccrs

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

def add_vertical_contours(w,lat,lon,
                          wmap_height=300, wmap_levels=[1,3],
                          annotate=True, xy0=[.01,-.05]):
    '''
    ARGS:
        w [lat,lon] : map of vertical motion
        lat,lon : degrees
        wmap_height : w map altitude
        annotate {True | False} : add annotation?
        xy0: annotate text added at xy0 and xy0 + [0,-.06]
    '''
    
    with warnings.catch_warnings():
        # ignore warning when there are no contours to plot:
        warnings.simplefilter('ignore')
        # downward motion in blueish
        plt.contour(lon, lat, -1*w, levels=wmap_levels, 
                    linestyles=['dashed','solid'],
                    colors=('aquamarine',))
        # upward motion in pink
        plt.contour(lon, lat, w, levels=wmap_levels, 
                    linestyles=['dashed','solid'],
                    colors=('pink',))
        
    plt.annotate(
        'vertical motion at %dm altitude'%wmap_height, 
        xy=xy0, 
        xycoords='axes fraction', 
        fontsize=8
        )
    plt.annotate(
        'dashed is %.1fm/s, solid is %.1fm/s'%(wmap_levels[0],wmap_levels[1]), 
        xy=[xy0[0],xy0[1]-.06], 
        xycoords='axes fraction', 
        fontsize=8
        )

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
            ff_lead = np.where(np.sum(ff<=0, axis=0)>0)[0][0]
            print("debug: fire front most western point is %.2f (=lon[%d])"%(lon[ff_lead],ff_lead))
            print("debug: transect goes from %.2f to %.2f"%(start[1],end[1]))
            
    
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
        add_vertical_contours(wmap,lat,lon,
                              wmap_height=wmap_height,
                              wmap_levels=[1,3],)
    
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
    print("DEBUG: slicex and slicez",slicex[0,:5],slicez[:15,0])
    
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
                       extent=[115.68, 116.15, -33.025,-32.79],
                       lats=None,lons=None, ff=None, sh=None, 
                       u10=None, v10=None, 
                       wmap=None, wmap_height=None,
                       topog=None, topog_contours=[50]):
    """
    Top down view of Waroona/Yarloop, adding fire front and heat flux and 10m winds
    ARGUMENTS:
        lats/lons are 1D arrays, required if adding other stuff
        ff [lats,lons] : firefront array
        sh [lats,lons] : sensible heat flux
        u10 [lats,lons] : 10m altitude longitudinal winds
        v10 [lats,lons] : 10m altitude latitudinal winds
        topog [lats,lons] : surface altitude
        topog_contours : list of values to contour topography at
    RETURNS:
        fig, ax, proj : proj is map projection used by tiff
    """
    # first create map from tiff file
    
    fig, ax, proj = plotting.map_tiff_qgis(
        fname="waroonaz_osm.tiff", 
        extent=extent,
        fig=fig,
        subplot_row_col_n=subplot_row_col_n,
        EPSG=4326, # SHOULD MATCH PLATECARREE
        )
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    
    if ff is not None:
        # add firefront
        cs_ff = plotting.map_fire(ff,lats,lons,transform=True)
    if sh is not None:
        # add hot spots for heat flux
        cs_sh, cb_sh = plotting.map_sensibleheat(sh,lats,lons, alpha=0.6)
    if u10 is not None:
        # winds, assume v10 is also not None        
        plt.streamplot(lons,lats,u10,v10, 
                       color='k', #transform=ccrs.PlateCarree(),
                       density=(0.6, 0.5))
        
        s10 = np.hypot(u10,v10)
        
        plotting.annotate_max_winds(s10, xytext=(0,1.02), s="10m wind max = %5.1f m/s")
    
    if wmap is not None:
        add_vertical_contours(wmap,lats,lons,
                              wmap_height=wmap_height,
                              wmap_levels=[1,3],
                              annotate=True,
                              xy0=[0.8,1.1])
    
    if topog is not None:
        ax.contour(lons,lats,topog,topog_contours,
                    colors='k', alpha=0.8, linewidths=1,
                    transform=ccrs.PlateCarree(),
                    )
    
    # set limits back to latlon limits
    ax.set_ylim(ylims[0],ylims[1])#, transform=ccrs.PlateCarree())  # outliers only
    ax.set_xlim(xlims[0],xlims[1])#, transform=ccrs.PlateCarree())
    
    return fig, ax, proj

def make_plots_emberstorm(model_run='waroona_run3', 
                          hours=None, 
                          extent=None,
                          wmap_height=300):
    """
    run emberstorm plotting method on model output read by my iris fio scripts
    """
    if extent is None:
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
            for transecti, transect in enumerate(_transects_[:5]):
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
    mr = 'waroona_run3'
    extent1 = [115.8, 116.1, -32.92,-32.82]
    transect1 = _transects_[0]
    
    hours=fio.model_outputs[mr]['filedates']
    testhours = [datetime(2016,1,6,13)]
    interesting_hours=[datetime(2016,1,6,x) for x in range(7,15)]

    if False:
        # This makes the combined 3 row plot with top down winds and 
        # transects of theta and wind
        make_plots_emberstorm(mr, hours=interesting_hours,
                              extent=extent1)
    
    if True:
        # newer plots showing 1: fire + town + winds (based on top panel in make_plots_emberstorm)
        # First of two emberstorms
        
        wmap_height = 300
        
        for extent,transect in zip([extent1,],[transect1,]):
                           
            dtimes=interesting_hours
            
            cubes = fio.read_model_run(mr, fdtime=dtimes, extent=extent, 
                                       add_topog=False, add_winds=True)
            w,=cubes.extract("upward_air_velocity")
            ctimes = utils.dates_from_iris(w)
            
            # extra vert map at ~ 300m altitude
            levh  = w.coord('level_height').points
            levhind = np.sum(levh<wmap_height)
            wmap = w[:,levhind]
            # read fire
            ff,sh,u10,v10 = fio.read_fire(model_run=mr,
                                          dtimes=ctimes, 
                                          extent=extent,
                                          sensibleheat=True,
                                          wind=True)
            
            lats = ff.coord('latitude').points
            lons = ff.coord('longitude').points
            for dti, dt in enumerate(ctimes):
                ffd = ff[dti].data.data
                shd = sh[dti].data.data
                u10d = u10[dti].data.data
                v10d = v10[dti].data.data
                wmapd = wmap[dti].data.data
                
                fig,ax,proj = topdown_emberstorm(extent=extent,
                                                 lats=lats,lons=lons,
                                                 ff=ffd, sh=shd, 
                                                 u10=u10d, v10=v10d,
                                                 wmap=wmapd,
                                                 wmap_height=wmap_height)
                
                ## Add dashed line to show where transect will be
                start,end =transect
                plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
                         linewidth=2, transform=ccrs.PlateCarree())
                
                ## Plot title
                LT = dt + timedelta(hours=8)
                plt.title(LT.strftime('%b %d, %H00(local)'))
                fio.save_fig(mr,_sn_,dt,subdir='topdown',plt=plt)
