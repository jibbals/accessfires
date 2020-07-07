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
import matplotlib.patheffects as PathEffects
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

# two extents worth examining:
_emberstorm_extents_ = [
        # first day around 10pm LT
        [115.8, 116.1, -32.92,-32.82], 
        # second day around 7pm LT
        [115.72, 116.03, -33.03, -32.839],
        ]
# lat lon pairs for start,end of transects
_emberstorm_transects_ = [
        [[-32.87,115.87],[-32.87,116.02]],
        [[-32.91,115.8],[-32.91,115.95]]
        ]

###
## METHODS
###

def add_vertical_contours(w,lat,lon,
                          wmap_height=300, wmap_levels=[1,3],
                          annotate=True, xy0=[.7,-.02], xy1=None):
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
    if xy1 is None:
        xy1=[xy0[0], xy0[1]-.03]
    plt.annotate(
        'vertical motion at %dm altitude'%wmap_height, 
        xy=xy0, 
        xycoords='axes fraction', 
        fontsize=8
        )
    plt.annotate(
        'dashed is %.1fm/s, solid is %.1fm/s'%(wmap_levels[0],wmap_levels[1]), 
        xy=xy1, 
        xycoords='axes fraction', 
        fontsize=8
        )

def transect_emberstorm(u,v,w,z,
                        lats,lons,
                        transect,
                        topog,
                        ztop=700,
                        ff=None,
                        theta=None,
                        theta_contourargs={},
                        ):
    """
    Plot transect showing east-west-north-south streamplot
    overlaid on horizontal wind contourf
    with potential temperature contours
    ARGUMENTS:
        u,v,w,z: arrays [lev,lats,lons] of East wind, N wind, Z wind, and level altitude
        lats,lons: degrees dimensions
        transect: [[lat0, lon0], [lat1,lon1]] transect (lat0 == lat1)
        topog: [lats,lons] surface altitude array
        ztop: top altitude to look at, defualt 800m
        ff: [lats,lons] firefront array (optional)
        theta: potential temperature (if contour is desired)
        theta_contourargs: dict(contour args for pot temp)
            defaults: levels=[300,305,310]
        
    """
    # First we subset all the arrays so to be below the z limit
    zmin = np.min(z,axis=(1,2)) # lowest altitude on each model level
    ztopi = np.argmax(ztop<zmin)+1 # highest index where ztop is less than model level altitude
    
    u,v,w,z = [u[:ztopi],v[:ztopi],w[:ztopi],z[:ztopi]]
    if theta is not None:
        theta=theta[:ztopi]
    start,end = transect
    # interpolation points
    npoints=utils.number_of_interp_points(lats,lons,start,end)
    # horizontal wind speed m/s
    s = np.hypot(u,v)
    # contourf of horizontal wind speeds
    #print("DEBUG:", s.shape, z.shape, lats.shape, lons.shape, start, end)
    _, slicex, slicez = plotting.transect_s(s,z, 
                                            lats,lons, 
                                            start, end,
                                            ztop=ztop,
                                            title="",
                                            lines=None, npoints=npoints,
                                            topog=topog, ff=ff)
    
    # east-west and vertical winds on transect
    sliceu,_,_ = utils.transect(u,lats,lons,start,end,nx=npoints, z_th=z)
    slicew,_,_ = utils.transect(w,lats,lons,start,end,nx=npoints, z_th=z)
    
    # Streamplot
    plotting.streamplot_regridded(slicex,slicez,sliceu,slicew,
                                  density=(1,1), 
                                  color='darkslategrey',
                                  zorder=1)
    plt.xlim(np.min(slicex),np.max(slicex))
    plt.ylim(np.min(slicez),ztop)
    
    
    ## Theta contours
    if theta is not None:
        
        sliceth,_,_ = utils.transect(theta,lats,lons,start,end,nx=npoints,z_th=z)
        # set defaults for theta contour plot
        if 'levels' not in theta_contourargs:
            theta_contourargs['levels'] = [295,300,305,310]
        if 'cmap' not in theta_contourargs:
            theta_contourargs['cmap'] = 'YlOrRd'#['grey','yellow','orange','red']
        if 'alpha' not in theta_contourargs:
            theta_contourargs['alpha'] = 0.9
        if 'linestyles' not in theta_contourargs:
            theta_contourargs['linestyles'] = 'dashed'
        if 'linewidths' not in theta_contourargs:
            theta_contourargs['linewidths'] = 0.9
        if 'extend' not in theta_contourargs:
            theta_contourargs['extend'] = 'both'
        
        # add faint lines for clarity
        contours = plt.contour(slicex,slicez,sliceth, **theta_contourargs)
        contours.set_clim(theta_contourargs['levels'][0], 
                          theta_contourargs['levels'][-1])
        plt.clabel(contours, inline=True, fontsize=10)
        
        
    
    

def emberstorm(theta, u, v, w, z, topog,
               lat, lon,
               ff=None,
               u10=None, v10=None,
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
    start,end=transect
    ## figure
    #fig, axes = plt.subplots(3,1)
    #plt.sca(axes[0])
    ## First plot, topography
    if u10 is None:
        u10=u[1]
    if v10 is None:
        v10=v[10]
    fig,ax0 = topdown_emberstorm(subplot_row_col_n=[3,1,1], extent=extent, lats=lat, lons=lon, ff=ff, sh=None, u10=u10, v10=v10, wmap=wmap, wmap_height=wmap_height, )
    ## THIS WAS ORIGINAL PLOTTING TOPDOWN VIEW
    #cs,cb = plotting.map_topography(extent,topog,lat,lon,title="Overview")
    #
    #start, end = transect
    #
    ## Add fire front contour
    #ff_lead_frac=None
    #if ff is not None:
    #    
    #    plotting.map_fire(ff,lat,lon)
    #    # west most burnt bit is fire front lead
    #    # grab first index from left to right along burnt longitudes
    #    if np.sum(ff<=0) > 0:
    #        #print("debug:", ff.shape, len(lat), len(lon), np.sum(ff<=0,axis=0).shape)
    #        ff_lead = np.where(np.sum(ff<=0, axis=0)>0)[0][0]
    #        print("debug: fire front most western point is %.2f (=lon[%d])"%(lon[ff_lead],ff_lead))
    #        print("debug: transect goes from %.2f to %.2f"%(start[1],end[1]))
    #        
    #
    ## start to end x=[lon0,lon1], y=[lat0, lat1]
    #plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
    #         linewidth=2,) 
    #         #marker='>', markersize=7, markerfacecolor='white')
    #
    ## Add markers where other transects will be 
    #if shadows is not None:
    #    for sstart,send in shadows:
    #        plt.plot([sstart[1],send[1]],[sstart[0],send[0], ], 
    #                 color='b', alpha=0.4, linestyle='None',
    #                 marker='>', markersize=2)
    #
    ## add nearby towns
    #plotting.map_add_locations(['waroona'], text=[''], color='grey',
    #                           marker='o', markersize=4)
    #plotting.map_add_locations(['fire_waroona'], text=[''], 
    #                           marker='*', color='r')
    #
    ## show winds
    #plt.streamplot(lon,lat,u[0],v[0], color='k',
    #               density=(.8,.5),
    #               )#alpha=0.7)
    #
    ## Finally add some faint pink or green hatching based on vertical wind motion
    #if wmap is not None:
    #    add_vertical_contours(wmap,lat,lon,
    #                          wmap_height=wmap_height,
    #                          wmap_levels=[1,3],)
    #
    ## cut down to desired extent
    #plt.ylim(extent[2:]); plt.xlim(extent[0:2])
    #
    ### Turn off the tick values
    #plt.xticks([]); plt.yticks([])
    
    ## Subplot 2, transect of potential temp
    # only looking up to 1km
    # horizontal points
    npoints = utils.number_of_interp_points(lat,lon,start,end)
    ax1=plt.subplot(3,1,2)
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
        ax1.annotate('', xy=(xslice[0,0],wmap_height), 
            xytext=(xslice[0,5], wmap_height),
            arrowprops=dict(facecolor='grey', arrowstyle='wedge,tail_width=0.5', alpha=0.5),
            fontsize=8)
    
    ## Next plot
    ax2=plt.subplot(3,1,3)
    transect_emberstorm(u,v,w,z,lat,lon,transect=[start,end],topog=topog,ztop=ztop,ff=ff)
    #_,slicex,slicez = plotting.transect_w(w,z,lat,lon,start,end,
    #                                      npoints=npoints,topog=topog, ff=ff,
    #                                      ztop=ztop)
    ## add east-west-vertical winds
    #uslice,_,_ = utils.transect(u,lat,lon,start,end,nx=npoints)
    ##vslice = utils.cross_section(v,lat,lon,start,end,npoints=npoints)
    #wslice,_,_ = utils.transect(w,lat,lon,start,end,nx=npoints)
    #
    ### Streamplot of east-west and vertical winds
    #
    #plotting.streamplot_regridded(slicex,slicez,uslice,wslice,
    #                              density=(.8,.5))
    
    return fig, [ax0,ax1,ax2]

def topdown_emberstorm(fig=None, subplot_row_col_n=None, 
                       extent=[115.68, 116.15, -33.025,-32.79],
                       lats=None,lons=None, ff=None, sh=None, 
                       u10=None, v10=None, 
                       wmap=None, wmap_height=None,
                       topog=None):
    """
    Top down view of Waroona/Yarloop, adding fire front and heat flux and 10m winds
    ARGUMENTS:
        lats/lons are 1D arrays, required if adding other stuff
        ff [lats,lons] : firefront array
        sh [lats,lons] : sensible heat flux
        u10 [lats,lons] : 10m altitude longitudinal winds
        v10 [lats,lons] : 10m altitude latitudinal winds
        topog [lats,lons] : surface altitude - can use this instead of tiff
        topog_contours : list of values to contour topography at - default 50m
    RETURNS:
        fig, ax
    """
    if fig is None:
        xsize = 12
        ysize = 12
        if extent is not None:
            # try to guess a good size for aspect ratio
            width = extent[1]-extent[0]
            height = extent[3]-extent[2]
            if width > (1.5*height):
                xsize=16
            if width > (2*height):
                xsize=20
                ysize=10
            if width < (0.75 * height):
                ysize=16
        fig=plt.figure(figsize=(xsize,ysize))
    # first create map from tiff file unless topography passed in
    if topog is None:
        fig, ax = plotting.map_tiff_qgis(
            fname="waroonaz_osm.tiff", 
            extent=extent,
            fig=fig,
            subplot_row_col_n=subplot_row_col_n,
            show_grid=True,
            aspect='equal',
            )
    else:
        plotting.map_topography(extent,topog,lats,lons,
                                cbar=False,title="")
        ax=plt.gca()
        ax.set_aspect('equal')
        ax.grid(color='gray',alpha=0.4) # add gridlines
        
        ## Add waroona, hamel, yarloop if possible
        for txt in ['Waroona','Hamel','Yarloop']:
            ax.annotate(txt, xy=np.array(plotting._latlons_[str.lower(txt)])[::-1],
                        xycoords="data", # lat lon xy as coords are platecarree
                        fontsize=12, ha="center",
                        color='k',
                        path_effects=[PathEffects.withStroke(linewidth=2,
                                                             foreground="w")])
            
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    
    if ff is not None:
        # add firefront
        cs_ff = plotting.map_fire(ff,lats,lons)
    if sh is not None:
        # add hot spots for heat flux
        cs_sh, cb_sh = plotting.map_sensibleheat(sh,lats,lons,alpha=0.6)
    if u10 is not None:
        # winds, assume v10 is also not None        
        plt.streamplot(lons,lats,u10,v10, 
                       color='k',
                       density=(0.6, 0.5))
        
        s10 = np.hypot(u10,v10)
        
        plotting.annotate_max_winds(s10, s="10m wind max = %5.1f m/s")
    
    if wmap is not None:
        add_vertical_contours(wmap,lats,lons,
                              wmap_height=wmap_height,
                              wmap_levels=[1,3],
                              annotate=True,
                              xy0=[0.6,-0.02])
    
    # set limits back to latlon limits
    ax.set_ylim(ylims[0],ylims[1])
    ax.set_xlim(xlims[0],xlims[1])
    
    return fig, ax

def explore_emberstorm(model_run='waroona_run3', 
                       hours=None, 
                       extent=None,
                       wmap_height=300,
                       ztop=800):
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
                #utcstamp = dtime.strftime("%b %d %H:%M (UTC)")
                ltstamp = (dtime+timedelta(hours=8)).strftime("%b %d %H:%M (UTC)")
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
                           ztop=ztop,
                           wmap=w[levhind].data, wmap_height=wmap_height,
                           transect=transect, shadows=shadows)
                
                # Save figure into folder with numeric identifier
                stitle="Emberstorm %s"%ltstamp
                plt.suptitle(stitle)
                distance=utils.distance_between_points(transect[0],transect[1])
                plt.xlabel("%.3fE - %.3fE = %.1fkm (at %.3fS)"%(transect[0][1], 
                           transect[1][1], distance/1e3, -1*transect[0][0]))
                fio.save_fig(model_run=model_run, script_name=_sn_, pname=dtime, 
                             plt=plt, subdir='transect%d'%transecti)
            
def zoomed_emberstorm_plots(hours=[0],
                            first=True,second=False,
                            topdown=True, transect=True,
                            topography=False,
                            wmap_height=300,):
    """
    Create zoomed in pictures showing top down winds and zmotion, with 
    either topography or open street maps underlay
    Arguments:
        hours: list of which hours to plot [0,...,23]
        first: True if plotting first emberstorm event
        second: True if plotting second emberstorm event
        topdown: True if topdown plot is desired
        transect: True if transect plot is desired
        topography: True if topography is desired instead of OSM
        wmap_height: how high to show vmotion contours, default 300m
    """
    dtimes=fio.model_outputs['waroona_run3']['filedates'][np.array(hours)]
    # extents: list of 4 element lists of E,W,S,N boundaries to look at
    # transects: list of [[lat,lon],[lat1,lon1]], transects to draw
    #    one for each extent
    extents=[]
    transects=[]
    subdirs=[]
    if first:
        extents.append(_emberstorm_extents_[0])
        transects.append(_emberstorm_transects_[0])
        subdirs=["first"]
    if second:
        # append second emberstorm scene
        extents.append(_emberstorm_extents_[1])
        transects.append(_emberstorm_transects_[1])
        subdirs.append("second")
    
    for extent,transect,subdir in zip(extents,transects,subdirs):
        
        cubes = fio.read_model_run(mr, fdtime=dtimes, extent=extent, 
                                   add_topog=True, add_winds=True,
                                   add_z=True, add_theta=True)
        u,v,w,z = cubes.extract(["u","v","upward_air_velocity","z_th"])
        theta, = cubes.extract("potential_temperature")
        topog=cubes.extract("surface_altitude")[0].data
        topogd = topog if topography else None
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
        zd = z.data.data
        for dti, dt in enumerate(ctimes):
            ffd = ff[dti].data.data
            LT = dt + timedelta(hours=8)
            
            if topdown:
                shd = sh[dti].data.data
                u10d = u10[dti].data.data
                v10d = v10[dti].data.data
                wmapd = wmap[dti].data.data
                
                fig,ax = topdown_emberstorm(extent=extent,
                                            lats=lats,lons=lons,
                                            ff=ffd, sh=shd, 
                                            u10=u10d, v10=v10d,
                                            topog=topogd,
                                            wmap=wmapd,
                                            wmap_height=wmap_height)
                
                ## Add dashed line to show where transect will be
                if transect is not None:
                    start,end =transect
                    plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
                             linewidth=2, alpha=0.7)
                
                ## Plot title
                plt.title(LT.strftime('%b %d, %H%M(local)'))
                fio.save_fig(mr,_sn_,dt,subdir=subdir+'/topdown',plt=plt)
            
            ## New plot for transect goes here
            if transect:
                transect_emberstorm(u[dti].data.data,
                                    v[dti].data.data,
                                    w[dti].data.data,
                                    zd,
                                    lats,lons,
                                    transect,
                                    topog=topog,
                                    ff=ffd,
                                    theta=theta[dti].data.data)
                
                plt.title(LT.strftime('Transect %b %d, %H%M(local)'))
                fio.save_fig(mr,_sn_,dt,subdir=subdir+'/transect',plt=plt)

if __name__ == '__main__':
    plotting.init_plots()
    mr = 'waroona_run3'
    extent1,extent2 = _emberstorm_extents_
    
    transect1 = _transects_[0]
    
    hours=fio.model_outputs[mr]['filedates']
    testhours = [datetime(2016,1,6,7)]
    interesting_hours=[datetime(2016,1,6,x) for x in range(7,15)]

    dtimes=interesting_hours[-5:] 

    if True:
        # This makes the combined 3 row plot with top down winds and 
        # transects of theta and wind
        explore_emberstorm(mr, hours=interesting_hours,
                           extent=extent1, ztop=3000)
    
    if True:
        # newer plots showing 1: fire + town + winds (based on top panel in make_plots_emberstorm)
        # First of two emberstorms
        zoomed_emberstorm_plots(
                first=False,
                second=True,
                topdown=True,
                transect=True,
                topography=True,
                hours=np.arange(36,39)
                )
        
