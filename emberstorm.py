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
from matplotlib.ticker import FormatStrFormatter
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
yspread = np.array([0, .02, -.02])
_transects_ = [ [[y,x0],[y,x1]] for y in y0+yspread ]


#-32.84, 115.93  # suburb centre: -32.8430, 115.8526
#-32.8725, 115.92

_emberstorm_centres_ = {
    'waroona_run3':{
        'first':{
            # where to centre the transect and UTC time
            'latlontime':[
                # start midway between hamel and waroona
                [-32.855,115.945,datetime(2016,1,6,3,30)], 
                # first point of interest soon after spot ignited
                [-32.875,115.96, datetime(2016,1,6,12,10)],
                [-32.88,115.96,datetime(2016,1,6,12,40)],
                [-32.865,115.925,datetime(2016,1,6,13,10)],
                # just before downslope run move to southern waroona
                [-32.86,115.92,datetime(2016,1,6,13,40)],
                [-32.86,115.92,datetime(2016,1,6,14,10)],
                # follow the run westwards a bit
                [-32.857,115.89,datetime(2016,1,6,14,40)],
                ],
            # W,E,S,N to look at
            'extent':[115.8, 116.1, -32.92,-32.82],
            # hours to look at
            'hours':np.arange(15,24),
            },
        'second':{
            ## NOW LOOKING TOWARDS SECOND OCCURRENCE:
            'latlontime':[
                # 7PM local time look at just north of yarloop
                [-32.95, 115.90,datetime(2016,1,7,11)],
                ],
            'extent':[115.72, 116.03, -33.03, -32.839],
            'hours':np.arange(31,40), 
            },
        },
    }

def emberstorm_centres(model_run, key, times, dx=0.07):
    """
        return list of latlons, interpolated to match where pcb are 
        spotted in model_run at times given by input times
    """
    latlontimes = _emberstorm_centres_[model_run][key]['latlontime']
    lats = [lat for lat,_,_ in latlontimes]
    lons = [lon for _,lon,_ in latlontimes]
    estimes = [time for _,_,time in latlontimes]
    
    # X is hours since 2015
    # interpolate lats and lons onto new list of datetimes
    es_X = [(dt - datetime(2015,1,1)).total_seconds()/3600.0 for dt in estimes]
    full_X = [(dt - datetime(2015,1,1)).total_seconds()/3600.0 for dt in times]
    full_lats = np.interp(full_X,es_X, lats, left=lats[0], right=lats[-1])
    full_lons = np.interp(full_X,es_X, lons, left=lons[0], right=lons[-1])
    #full_latlons = [ (lat, lon) for lat,lon in zip(full_lats, full_lons) ]
    full_transects = [[[lat,lon-dx],[lat,lon+dx]] for lat,lon in zip(full_lats,full_lons)]
    return full_transects

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
                        sh=None,
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
                                            topog=topog, sh=sh)
    
    # east-west and vertical winds on transect
    sliceu,_,_ = utils.transect(u,lats,lons,start,end,nx=npoints, z_th=z)
    slicew,_,_ = utils.transect(w,lats,lons,start,end,nx=npoints, z_th=z)
    
    # Streamplot
    plotting.streamplot_regridded(slicex,slicez,sliceu,slicew,
                                  density=(1,1), 
                                  color='darkslategrey',
                                  zorder=1,
                                  #linewidth=np.hypot(sliceu,slicew), # too hard to see what's going on
                                  )
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
        if subplot_row_col_n is not None:
            prow,pcol,pnum=subplot_row_col_n
            ax = plt.subplot(prow,pcol,pnum)
        plotting.map_topography(extent,topog,lats,lons,
                                cbar=False,title="")
        ax=plt.gca()
        ax.set_aspect('equal')
        
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
        plt.annotate(s="max heat flux = %6.1e W/m2"%np.max(sh),
                     xy=[0,1.11],
                     xycoords='axes fraction', 
                     fontsize=8)
    if u10 is not None:
        # winds, assume v10 is also not None
        s10 = np.hypot(u10,v10)
        lw10 = utils.wind_speed_to_linewidth(s10, lwmax=5, speedmax=20)
        # higher density if using topography instead of OSM
        density=(0.6,0.5) if topog is None else (0.75,0.7)
        plt.streamplot(lons,lats,u10,v10, 
                       linewidth=lw10, 
                       color='k',
                       density=density,
                       )
        plotting.annotate_max_winds(s10, s="10m wind max = %5.1f m/s")
    
    if wmap is not None:
        add_vertical_contours(wmap,lats,lons,
                              wmap_height=wmap_height,
                              wmap_levels=[1,3],
                              annotate=True,
                              xy0=[0.7,1.08])
        
    # set limits back to latlon limits
    ax.set_ylim(ylims[0],ylims[1])
    ax.set_xlim(xlims[0],xlims[1])
    # 115.8, 116.1, -32.92,-32.82
    xticks=np.arange(115.8,116.11,0.05)
    plt.xticks(xticks,xticks)
    yticks=np.arange(-32.92,-32.805,0.03)
    plt.yticks(yticks,yticks)
    
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    # add gridlines
    #gl=ax.grid(color='gray',alpha=0.4) 
    #gl.xlabels_top = False
    #gl.xlabels_bottom = True
    #gl.ylabels_left = True
    #gl.ylabels_right = False
    
    return fig, ax

def explore_emberstorm(model_run='waroona_run3', 
                       hours=None, 
                       extent=None,
                       topography=False,
                       wmap_height=300,
                       ztop=800):
    """
    run emberstorm plotting method on model output read by my iris fio scripts
    ARGUMENTS:
        hours=[datetimes]
        extent=[W,E,S,N]
        topography=False, set true to use topog for topdown view
        wmap_height=300, what height for topdown vertical motion contours?
        ztop=800, how high to do transect?
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
        ff, sh, u10, v10 = fio.read_fire(model_run=model_run, dtimes=dtimes, 
                                     extent=extent,
                                     wind=True, sensibleheat=True)
        
        # pull out bits we want
        uvw = cubelist.extract(['u','v','upward_air_velocity'])
        # extra vert map at ~ 300m altitude
        
        levh  = uvw[2].coord('level_height').points
        levhind = np.sum(levh<wmap_height)
        
        z, = cubelist.extract(['z_th']) # z has no time dim
        # for each time slice pull out potential temp, winds
        for i,dtime in enumerate(dtimes):
            for transecti, transect in enumerate(_transects_):
                
                #utcstamp = dtime.strftime("%b %d %H:%M (UTC)")
                ltstamp = (dtime+timedelta(hours=8)).strftime("%b %d %H:%M (LT)")
                # winds
                u,v,w = uvw[0][i].data.data, uvw[1][i].data.data, uvw[2][i].data.data
                
                # fire
                ffi,shi=None,None
                if ff is not None:
                    ffi = ff[i].data 
                if sh is not None:
                    shi = sh[i].data
                #vertical motion at roughly 300m altitude
                wmap=w[levhind]
                topogd=topog.data if topography else None
                
                start,end=transect
                
                ## First plot, topography
                fig,ax0 = topdown_emberstorm(
                        subplot_row_col_n=[2,1,1], 
                        extent=extent, lats=lat, lons=lon, 
                        topog=topogd,
                        ff=ffi, sh=shi, 
                        u10=u10[i].data, v10=v10[i].data,
                        wmap=wmap, wmap_height=wmap_height, 
                        )
                
                ## Add transect line
                # start to end x=[lon0,lon1], y=[lat0, lat1]
                plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
                         linewidth=2, 
                         #marker='>', markersize=7, markerfacecolor='white'
                         )
                
                ## Subplot 2, transect of potential temp
                # how many horizontal points to interpolate to
                #npoints = utils.number_of_interp_points(lat,lon,start,end)
                #ax1=plt.subplot(3,1,2)
                #trets = plotting.transect_theta(theta[i].data, z.data, lat, lon, start, end, npoints=npoints,
                #                                topog=topog.data, ff=ffi, ztop=ztop,
                #                                contours=np.arange(290,320.1,0.5),
                #                                lines=None, #np.arange(290,321,2), 
                #                                linestyles='dashed')
                ## add faint lines for clarity
                #thetaslice,xslice,zslice=trets
                #plt.contour(xslice,zslice,thetaslice,np.arange(290,320.1,1),colors='k',
                #            alpha=0.5, linestyles='dashed', linewidths=0.5)
                # 
                ## Finally show winds on transect
                ax2=plt.subplot(2,1,2)
                
                transect_emberstorm(u,v,w,z.data,lat,lon,transect=[start,end],topog=topog.data,ztop=ztop,sh=shi)
                
                # Save figure into folder with numeric identifier
                stitle="Emberstorm %s"%ltstamp
                plt.suptitle(stitle)
                distance=utils.distance_between_points(transect[0],transect[1])
                plt.xlabel("%.3fE - %.3fE = %.1fkm (at %.3fS)"%(transect[0][1], 
                           transect[1][1], distance/1e3, -1*transect[0][0]))
                fio.save_fig(model_run=model_run, script_name=_sn_, pname=dtime, 
                             plt=plt, subdir='transect%d'%transecti)
            
def zoomed_emberstorm_plots(hours=None,
                            first=True,second=False,
                            topography=False,
                            extent=None,
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
    for i in range(2):
        if (i==0) and not first:
            continue
        if (i==1) and not second:
            continue
        key=['first','second'][i]
        
        if hours is None:
            hours = _emberstorm_centres_['waroona_run3'][key]['hours']
        if extent is None:
            extent = _emberstorm_centres_['waroona_run3'][key]['extent']
        
        dtimes=fio.model_outputs['waroona_run3']['filedates'][np.array(hours)]
        
        cubes = fio.read_model_run(mr, fdtime=dtimes, extent=extent, 
                                   add_topog=True, add_winds=True,
                                   add_z=True, add_theta=True)
        u,v,w,z = cubes.extract(["u","v","upward_air_velocity","z_th"])
        theta, = cubes.extract("potential_temperature")
        topog=cubes.extract("surface_altitude")[0].data
        topogd = topog if topography else None
        ctimes = utils.dates_from_iris(w)
        
        # transects: list of [[lat,lon],[lat1,lon1]], transects to draw
        transects=emberstorm_centres('waroona_run3',key,ctimes)
        
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
            transect=transects[dti]
            shd = sh[dti].data.data
            LT = dt + timedelta(hours=8)
            
            ffd = ff[dti].data.data
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
            fio.save_fig(mr,_sn_,dt,subdir=key+'/topdown',plt=plt)
            
            ## New plot for transect goes here
            transect_emberstorm(u[dti].data.data,
                                v[dti].data.data,
                                w[dti].data.data,
                                zd,
                                lats,lons,
                                transect,
                                topog=topog,
                                sh=shd,
                                theta=theta[dti].data.data)
            
            plt.title(LT.strftime('Transect %b %d, %H%M(local)'))
            fio.save_fig(mr,_sn_,dt,subdir=key+'/transect',plt=plt)

if __name__ == '__main__':
    plotting.init_plots()
    mr = 'waroona_run3'
    #extent1,extent2 = _emberstorm_extents_
    extent1 = _emberstorm_centres_['waroona_run3']['first']['extent']
    extent2 = _emberstorm_centres_['waroona_run3']['second']['extent']
    #transect1 = _transects_[0]
    
    hours=fio.model_outputs[mr]['filedates']
    testhours = [datetime(2016,1,6,7)]
    interesting_hours=hours[12:24]

    dtimes=interesting_hours[-7:] 

    if False:
        # This makes the combined 3 row plot with top down winds and 
        # transects of theta and wind
        # Let's do half with topography and half with OSM
        #explore_emberstorm(mr, hours=hours[7:13],
        #                   topography=False,
        #                   extent=extent1, ztop=3000)
        explore_emberstorm(mr, hours=hours[18:24],
                           topography=True,
                           extent=extent1, ztop=3000)
    
    if True:
        # newer plots showing 1: fire + town + winds (based on top panel in make_plots_emberstorm)
        # First of two emberstorms
        zoomed_emberstorm_plots(
                first=True,
                second=True,
                topography=True,
                #hours=np.arange(15,24)
                )
        
