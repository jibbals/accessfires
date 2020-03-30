# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 11:09:54 2019
    Zoom in on pyrocb
    Plot:
        311: top-down view with transect[s] plotted as lines
        312: vert motion transect along emberstorm
        313: vert motion transect along emberstormx (roughly perp to emberstorm)

@author: jgreensl
"""

import matplotlib
matplotlib.use('Agg',warn=False)

# plotting stuff
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.ticker as tick
import cartopy.crs as ccrs

import numpy as np
import pandas as pd
from datetime import datetime,timedelta
import iris # file reading and constraints etc
import warnings

# local modules
from utilities import plotting, utils, fio, constants

###
## GLOBALS
###
_sn_ = 'pyrocb'

__cloud_thresh__ = constants.cloud_threshold


## List of PyroCB latlon, and time 
__PCB_occurrences__ = {
    'waroona_run3':{
        'latlon':[(-32.9,116.05), (-32.9,116.04),
                  (-32.9, 116.00), (-32.9,115.96),
                  (-32.88, 116.01), (-32.88,115.94)],
        'time':[datetime(2016,1,6,6), datetime(2016,1,6,7), 
                datetime(2016,1,6,8), datetime(2016,1,6,10),
                datetime(2016,1,6,11),datetime(2016,1,6,14)]
        },
    'waroona_run2':{
        'latlon':[(-32.9,116.08),(-32.91,116.07),
                  (-32.9,116.02),(-32.89,116.00),
                  (-32.9,115.99),(-32.90,115.97),
                  (-33.0,115.88), (-32.88,115.94)],
        'time':[datetime(2016,1,6,5), datetime(2016,1,6,6), 
                datetime(2016,1,6,7), datetime(2016,1,6,8),
                datetime(2016,1,6,9), datetime(2016,1,6,10),
                datetime(2016,1,6,11),datetime(2016,1,6,14)]
        },
    'waroona_run1':{
        'latlon':[(-32.88,116.14),(-32.88,116.0),],
        'time':[datetime(2016,1,5,15), datetime(2016,1,6,14),]
        },
    'waroona_old':{
        'latlon':[(-32.88,116.14), (-32.89,116.14),        
                  (-32.9,116.11),],
        'time':[datetime(2016,1,5,21), datetime(2016,1,6,0), 
                datetime(2016,1,6,8),]
        },
    'sirivan_run2_hr':{
        'latlon':[(-32.1,149.725),(-32.1,149.875),
                  (-32.04,149.875),(-32.05,150.0255),
                  (-31.9,150)],
        'time':[datetime(2017,2,12,1), datetime(2017,2,12,3), 
                datetime(2017,2,12,6), datetime(2017,2,12,7),
                datetime(2017,2,12,18)]
        },
    'sirivan_run2':{
        'latlon':[(-32.1,149.725),(-32.1,149.875),
                  (-32.04,149.875),(-32.05,150.0255),
                  (-31.9,150)],
        'time':[datetime(2017,2,12,1), datetime(2017,2,12,3), 
                datetime(2017,2,12,6), datetime(2017,2,12,7),
                datetime(2017,2,12,18)]
        },
    'sirivan_run1':{
        'latlon':[(-32.1,149.725),(-32.1,149.875),
                  (-32.04,149.875),(-32.05,150.0255),
                  (-31.9,150)],
        'time':[datetime(2017,2,12,1), datetime(2017,2,12,3), 
                datetime(2017,2,12,6), datetime(2017,2,12,7),
                datetime(2017,2,12,18)]
        },
    }

def pcb_occurrences(model_run, times):
    """
        return list of latlons, interpolated to match where pcb are 
        spotted in model_run at times given by input times
    """
    latlons = __PCB_occurrences__[model_run]['latlon']
    lats = [lat for lat,_ in latlons]
    lons = [lon for _,lon in latlons]
    pcbtimes = __PCB_occurrences__[model_run]['time']
    
    # X is hours since 2015
    # interpolate lats and lons onto new list of datetimes
    pcb_X = [(dt - datetime(2015,1,1)).total_seconds()/3600.0 for dt in pcbtimes]
    full_X = [(dt - datetime(2015,1,1)).total_seconds()/3600.0 for dt in times]
    full_lats = np.interp(full_X,pcb_X, lats, left=lats[0], right=lats[-1])
    full_lons = np.interp(full_X,pcb_X, lons, left=lons[0], right=lons[-1])
    full_latlons = [ (lat, lon) for lat,lon in zip(full_lats, full_lons) ]
    
    return full_latlons

def transect_plus_quiver(w,u,qc,topog,zth,lat,lon,start,end,ztop,contours,
                         npoints=100,nquivers=13):
    '''
    '''
    wslice, xslice, zslice = plotting.transect_w(w, zth,
                                                 lat, lon, start, end,
                                                 npoints=npoints, title='',
                                                 topog=topog, ztop=ztop,
                                                 contours=contours,
                                                 lines=None, colorbar=True)
    ## add cloud outlines
    ## Add contour where clouds occur
    qcslice = utils.cross_section(qc, lat, lon, start, end, npoints=npoints)
    with warnings.catch_warnings():
        # ignore warning when there are no clouds:
        warnings.simplefilter('ignore')
        plt.contour(xslice,zslice,qcslice,np.array([__cloud_thresh__]),colors='k')

    
    ## Add quivers to transect
    # get longitudinal wind speed slice along transect
    uslice = utils.cross_section(u,lat,lon,start,end,npoints=npoints)
    
    # lets cut away the upper levels before making our quiver
    ztopirows, ztopicols = np.where(zslice < ztop) # index rows and columns
    ztopi = np.max(ztopirows) # highest index with height below ztop
    
    # quiver east-west and vertical winds
    quiverscale=160
    plotting.map_quiver(uslice[:ztopi+4,:], wslice[:ztopi+4,:], 
                        zslice[:ztopi+4,:], xslice[:ztopi+4,:], 
                        nquivers=nquivers, scale=quiverscale,
                        alpha=0.5, pivot='middle')


def map_with_transect(data,lat,lon, transect,
                      ff=None, 
                      color='m', linestyle='--', linewidth=2,
                      extralines=[], extracolors=[],
                      **map_contourf_args):
    """
    plot contourf of whatever, plus firefront (optional).
    Add transect lines on top of contourf.
    """
    extent=[lon[0],lon[-1],lat[0],lat[-1]]
    start,end = transect
    
    # map contour using input data and arguments
    cs, cb = plotting.map_contourf(extent, data, lat, lon,
                                   **map_contourf_args)
    
    # fn for adding little arrow to map
    def myarrow(yx0, yx1):
        y0,x0=yx0
        y1,x1=yx1
        dx,dy = (x1-x0)/15.0,(y1-y0)/15.0
        plt.arrow(x0-1.4*dx,y0-1.4*dy,dx,dy,width=0.0015)
    
    # start to end x=[lon0,lon1], y=[lat0, lat1]
    plt.plot([start[1],end[1]],[start[0],end[0], ], linestyle+color, 
             linewidth=linewidth)
    myarrow(start,end)
    
    # any extra lines do them too
    for [startx,endx],xcolor in zip(extralines, extracolors):
        plt.plot([startx[1],endx[1]],[startx[0],endx[0], ], '--',
                 color=xcolor,
                 linewidth=2)
        myarrow(startx,endx)
    
    # Add fire outline
    if ff is not None:
        plotting.map_fire(ff,lat,lon)
    
    return cs,cb

def left_right_slice(qc, u, w, z,
                     topog, lat, lon, 
                     transect,
                     ztop=4000,
                     cloud_threshold=constants.cloud_threshold,
                     npoints=100, nquivers=15
                     ):
    '''
    cloud content contourf, with u,w wind quiver overlaid
    
    INPUTS:
        qc: cloud liquid + ice g/kg [z, lat, lon]
        u,w : horizontal, vertical winds m/s [z,lat,lon]
        z: model level heights m [z, lat, lon]
        topog: surface altitude m [lat,lon]
        transect: [[lat0,lon0],[lat0,lon1]]
        ztop: height of z axis in transect
        cloud_threshold: g/kg threshold of liquid + ice in air
        npoints: how many points to interpolate to on the transect
    '''
    start,end = transect
       
    ## Plot vert motion transect
    qcslice, xslice, zslice = plotting.transect_qc(qc, z, lat, lon, start, end,
                                                   npoints=npoints, title='',
                                                   topog=topog, ztop=ztop,
                                                   colorbar=False)
    ## add wind quivers
    # need transect u and w
    uslice = utils.cross_section(u, lat, lon, start, end, npoints=npoints)
    wslice = utils.cross_section(w, lat, lon, start, end, npoints=npoints)
    
    # quiver without crazy density
    vsx,vsz = np.array(xslice.shape)//nquivers
    # also skip first few so that barbs don't extend the transect
    skip = (slice(vsx,None,vsx),slice(None,None,vsz))
    
    # meters per second to knots
    mpsk = 1.94384
    plt.barbs(xslice[skip],zslice[skip],uslice[skip]*mpsk,wslice[skip]*mpsk, 
              pivot='middle')
    #plt.quiver(xslice[skip],zslice[skip],uslice[skip],wslice[skip], 
    #           scale=quiverscale, pivot='Middle')
    
    plt.ylabel('height (m)')
    plt.xlabel('')
    
def top_down_vertical_motion(W,lat,lon,FF=None,Q=None,):
    """
    W[lat,lon]: vertical motion (m/s)
    FF, Q : firefront, cloud content(g/kg)
    """
    
    # Standard color scale
    wcmap=plotting._cmaps_['verticalvelocity']
    wnorm=col.SymLogNorm(0.25) # linear to +- 0.25, then log scale
    wcontours=np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0]))
    
    # Vertical motion on standard colorscale
    cs = plt.contourf(lon, lat, W, wcontours, cmap=wcmap, norm=wnorm)

    # add cloud hatching
    if Q is not None:
        plt.contourf(lon,lat,Q, [0,__cloud_thresh__,100], 
                     colors=['None']*3, hatches=[None,'/','//'],
                     extend='both',)
        
    if FF is not None:
        plotting.map_fire(FF,lat,lon)
        
    plt.xticks([],[])
    plt.yticks([],[])
        
    return cs    
    ## reduce vert gap between subplots
    #fig.subplots_adjust(hspace=0.1)
    ## add vert wind colourbar
    #cbar_ax = fig.add_axes([0.905, 0.4, 0.01, 0.2])# X Y Width Height
    #fig.colorbar(cs, cax=cbar_ax, format=ticker.ScalarFormatter(), pad=0)

def sample_showing_grid(model_run='waroona_run3', extentname=None, HSkip=None):
    """
    Show each hour the latlon grid and vertical motion at 4000,5000,6000 metres
    """
    all_hours = fio.model_outputs[model_run]['filedates']
    
    if extentname is None:
        extentname=model_run.split('_')[0]
    
    extent=plotting._extents_[extentname]
    
    for di,hour in enumerate(all_hours):
        ## Read hour of output
        try:
            cubes=fio.read_model_run(model_run,fdtime=hour,extent=extent,
                                     add_winds=True, HSkip=HSkip)
        except OSError as ose:
            print("WARNING: OSError: probably missing some data")
            print("       :", ose)
            continue
        
        W, = cubes.extract('upward_air_velocity')
        Q, = cubes.extract('qc')
        height = W.coord('level_height').points
        lat = W.coord('latitude').points
        lon = W.coord('longitude').points
        
        h1 = np.argmin(np.abs(height-4000)) # nearest to 4k
        h2 = np.argmin(np.abs(height-5000))
        h3 = np.argmin(np.abs(height-7000))
        FF, = fio.read_fire(model_run, dtimes=[hour], extent=extent, HSkip=HSkip)
        ## Horizontal map of vertical motion
        f, axes = plt.subplots(3,1,figsize=[12,16],subplot_kw={'projection': ccrs.PlateCarree()})
        for ax, hi in zip(axes,[h1,h2,h3]):
            plt.sca(ax)
            cs = top_down_vertical_motion(W[0,hi].data.data,
                                          FF=FF[0].data, Q=Q[0,hi].data.data,
                                          lat=lat,lon=lon,)
            plotting.map_add_locations_extent(extentname,hide_text=True)
            plotting.map_draw_gridlines(ax,)
            #plt.setp(ax.get_xticklabels(), rotation=45)
            plt.title('~%d m'%int(height[hi]))
        
        fio.save_fig(model_run,_sn_,hour,plt,subdir="sample_with_grid")
        

def pyrocb(w, u, qc, z, wmean, topog, lat, lon,
           transect1, transect2, transect3,
           ff=None,
           wmeantitle='Average vertical motion',
           extentname=None,
           ztop=15000,
           cloud_threshold=constants.cloud_threshold,
           ):
    '''
    311: mean vert windspeed between 500m and 1500m (?)
    312: vert motion transect (long transect)
    325: vert motion transect (short transect 1)
    326: vert motion transect (short transect 2)
    
    INPUTS:
        w: vertical motion m/s [z, lat, lon]
        qc: cloud liquid + ice g/kg [z, lat, lon]
        z: model level heights m [z, lat, lon]
        wmean: mean vertical motion for subplot 311 m/s [lat,lon]
        topog: surface altitude m [lat,lon]
        transectN: [[lat0,lon0],[lat1,lon1]] of transect
        ff: firefront array (-ve is burnt), optional
        ztop: height of z axis in transect
        cloud_threshold: g/kg threshold of liquid + ice in air
        ext is the plot extension { '.png' | '.eps' }
        dpi is plot quality (100 is default, higher is better for publication)
    '''
    ## Plotting setup
    # set font sizes etc
    plotting.init_plots()
    # get transects
    start,end = transect1
    startx1,endx1 = transect2
    startx2,endx2 = transect3
    colorx1='b'
    colorx2='teal'
    
       
    # figure setup
    f=plt.figure(figsize=[7,10])
    plt.subplot(3,1,1)

    ### Plot 311
    # top panel is wind speed surface values
    clevs_vertwind = np.union1d(np.union1d(2.0**np.arange(-2,6),
                                           -1*(2.0**np.arange(-2,6))),
                                np.array([0]))
    cs, _ = map_with_transect(wmean, lat, lon, transect1, ff=ff, 
                              extralines=[transect2,transect3],
                              extracolors=['b','teal'],
                              clevs = clevs_vertwind, cbar=False, 
                              cmap=plotting._cmaps_['verticalvelocity'],
                              norm=col.SymLogNorm(0.25), 
                              cbarform=tick.ScalarFormatter(), clabel='m/s')
    plt.title(wmeantitle)
    
    # add nearby towns
    if extentname is not None:
        plotting.map_add_locations_extent(extentname)
    
    ### transect plots
    ###
    ax2=plt.subplot(3,1,2)
    
    ## Plot vert motion transect
    #    wslice, xslice, zslice = plotting.transect_w(w, z, lat, lon, start, end,
    #                                                 npoints=100, title='',
    #                                                 topog=topog, ztop=ztop,
    #                                                 lines=None, colorbar=False)
    #    ## add cloud outlines
    #    ## Add contour where clouds occur
    #    qcslice = utils.cross_section(qc, lat, lon, start, end, npoints=100)
    #    with warnings.catch_warnings():
    #        # ignore warning when there are no clouds:
    #        warnings.simplefilter('ignore')
    #        plt.contour(xslice,zslice,qcslice,np.array([cloud_threshold]),colors='k')
    transect_plus_quiver(w,u,qc,topog,z,lat,lon,start,end,ztop,
                         contours=clevs_vertwind,
                         npoints=100,nquivers=13)
    
    plt.ylabel('height (m)')
    plt.xlabel('')
    for spine in [ax2.spines['bottom'], ax2.spines['top']]:
        spine.set_linestyle('dashed')
        spine.set_capstyle("butt")
        spine.set_color("m")
        spine.set_linewidth(2)
    
    ax31, ax32 = [plt.subplot(3,2,5), plt.subplot(3,2,6)]
    for ax3,startx,endx,colorx in zip([ax31, ax32],
                                      [startx1,startx2],
                                      [endx1,endx2],
                                      [colorx1,colorx2]):
        plt.sca(ax3)
    
        ## Plot vert motion transect
        wslicex,xslicex,zslicex = plotting.transect_w(w, z, lat, lon, 
                                                      startx, endx,
                                                      title='',
                                                      npoints=100,
                                                      topog=topog, 
                                                      colorbar=False,
                                                      ztop=ztop,
                                                      lines=None)
        ## add cloud outlines
        ## Add contour where clouds occur
        qcslicex = utils.cross_section(qc, lat, lon, startx, endx, 
                                       npoints=100)
        with warnings.catch_warnings():
            # ignore warning when there are no clouds:
            warnings.simplefilter('ignore')
            plt.contour(xslicex,zslicex,qcslicex,np.array([cloud_threshold]),colors='k')

        plt.xlabel('')
        for spine in [ax3.spines['bottom'], ax3.spines['top']]:
            spine.set_linestyle('dashed')
            spine.set_capstyle("butt")
            spine.set_color(colorx)
            spine.set_linewidth(2)
    
    # remove y ticks from last figure
    plt.yticks([],[])
        
    # add colour bar
    axes=[0.85, 0.20, 0.03, 0.6] 
    #f.subplots_adjust(top=0.95)
    #f.subplots_adjust(hspace=0.05)
    f.subplots_adjust(right=axes[0]-0.01)
    cbar_ax = f.add_axes(axes)
    cb = f.colorbar(cs, cax=cbar_ax, 
                    format=tick.ScalarFormatter())
    # -ve labelpad moves label closer to colourbar
    cb.set_label('m/s', labelpad=-3)

def moving_pyrocb(model_run='waroona_run3', hours = None,
                  ztop=14000, xlen=0.1, nquivers=12):
    """
    follow pyrocb with a transect showing vert motion and pot temp
    ARGUMENTS:
        model_run: name of model run to look at
        hours: list of hours to plot
        ztop: altitute ceiling (metres) for transects
        xlen: horizontal distance (degrees) shown by transects
    """
    
    extentname=model_run.split('_')[0]
    extent = plotting._extents_[extentname]
    clevs_vertwind = np.union1d(np.union1d(2.0**np.arange(-2,6),
                                           -1*(2.0**np.arange(-2,6))),
                                np.array([0]))
    ztop=14000
    
    if hours is None:
        hours = fio.model_outputs[model_run]['filedates']

    ## read um output over extent [t, lev, lat, lon]
    for hour in hours:
        cubes = fio.read_model_run(model_run, extent=extent, 
                                   fdtime=[hour],
                                   add_topog=True,
                                   add_winds=True)
                                   #add_z=True)
        w, = cubes.extract('upward_air_velocity')
        ffdtimes = utils.dates_from_iris(w) 
        u, = cubes.extract('u')
        qc, = cubes.extract('qc')
        topog0, = cubes.extract('surface_altitude')
        topog = topog0.data
        #z, = cubes.extract('z_th') 
        
        lat = w.coord('latitude').points
        lon = w.coord('longitude').points
        
        pcb_centres = pcb_occurrences(model_run,times=ffdtimes)
        # transects = [lat0,lon0,lat1,lon1]
        transects = [[plat,plon-xlen/2.0,plat,plon+xlen/2.0] for (plat,plon) in pcb_centres]
        
        ## fire front
        ff, = fio.read_fire(model_run=model_run, dtimes=ffdtimes, extent=extent, firefront=True)
        # add zth cube
        p, pmsl = cubes.extract(['air_pressure','air_pressure_at_sea_level'])
        Ta, = cubes.extract('air_temperature')
        nz,ny,nx = p[0].shape
        
        for i,transect in enumerate(transects):
            # take mean of vert motion between lvls 25-48 approx 500m - 1500m
            with warnings.catch_warnings():
                # ignore warning from collapsing non-contiguous dimension:
                warnings.simplefilter('ignore')
                wmean = w[i,25:48,:,:].collapsed('model_level_number', iris.analysis.MEAN)
            h0,h1 = wmean.coord('level_height').bounds[0]
            fire=None
            if ff is not None:
                fire=ff[i].data.data
            qci = qc[i].data.data
            ui = u[i].data.data
            wi = w[i].data.data
            
            ## calc zth
            # repeat surface pressure along z axis
            reppmsl = np.repeat(pmsl[i].data[np.newaxis,:,:],nz, axis=0)
            zth = -(287*300/9.8)*np.log(p[i].data/reppmsl)
            
            a,b,c,d = transect
            start=[a,b]
            end=[c,d]
            
            ## map showing transect
            plt.figure(figsize=[10,16])
            plt.subplot(3,1,1)
            cs, _ = map_with_transect(wmean.data, lat, lon, transect=[start,end], 
                                      ff=fire, color='k', linewidth=1,
                                      clevs = clevs_vertwind,
                                      cbar=False, 
                                      cmap=plotting._cmaps_['verticalvelocity'],
                                      norm=col.SymLogNorm(0.25), 
                                      cbarform=tick.ScalarFormatter(), 
                                      clabel='m/s')
            plotting.map_add_locations_extent(extentname,hide_text=True) 
            wmeantitle='Vertical motion mean (%3.0fm - %4.0fm)'%(h0,h1)
            plt.title(wmeantitle)
            
            ## Transect of vert motion
            plt.subplot(3,1,2)
            transect_plus_quiver(wi,ui,qci,topog,zth,lat,lon,start,end,ztop,
                                 contours=clevs_vertwind,
                                 npoints=100,nquivers=13)
            
            plt.ylabel('height (m)')
            plt.xlabel('')
            plt.title('vertical motion transect')
            
            ## plot potential temp
            plt.subplot(3,1,3)
            theta = utils.potential_temperature(p[i].data,Ta[i].data)
            plotting.transect_theta(theta,zth,lat,lon,start,end,
                                    npoints=100, topog=topog, title='',
                                    ztop=ztop)
            
            ## Plot title and vmotion colour bar
            plt.title('potential temperature transect')
            stitle = ffdtimes[i].strftime("Vertical motion %Y %b %d %H:%M (UTC)")
            plt.suptitle(stitle)
            
            #fig=plt.gcf()
            #axes=[0.33, 0.65, 0.33, 0.02]
            #f.subplots_adjust(right=axes[0]-0.01)
            #cbar_ax = fig.add_axes(axes)
            #cb = fig.colorbar(cs, cax=cbar_ax, 
            #                format=tick.ScalarFormatter())
            # -ve labelpad moves label closer to colourbar
            #cb.set_label('m/s', labelpad=-3)
            
            fio.save_fig(model_run,_sn_,ffdtimes[i],plt,subdir='moving')
        
    
def pyrocb_model_run(model_run='waroona_run1', dtime=datetime(2016,1,5,15)):
    """
    Try to show pyrocb using two figures:
        1: left to right transect showing winds and potential temp
        2: Three transects (forming somewhat of an asterix) of vertical winds
    Makes figures for single hour defined by dtime input argument
    """
    ### First use datetime and extentname to read correct outputs:
    extentname=model_run.split('_')[0]
    extent = plotting._extents_[extentname]
    
    ## read um output over extent [t, lev, lat, lon]
    cubes = fio.read_model_run(model_run, fdtime=[dtime], extent=extent,
                               add_z=True, add_winds=True, add_topog=True)
    
    w, = cubes.extract('upward_air_velocity')
    u, = cubes.extract('u')
    qc, = cubes.extract('qc')
    topog, = cubes.extract('surface_altitude')
    z, = cubes.extract('z_th') 
    
    lat = w.coord('latitude').points
    lon = w.coord('longitude').points

    ## fire front
    ffdtimes = utils.dates_from_iris(w)
    ff, = fio.read_fire(model_run=model_run, dtimes=ffdtimes, extent=extent, firefront=True)
    
    lat_pcb, lon_pcb = pcb_occurrences(model_run,ffdtimes)[0]
    multiplyer=1
    if 'sirivan' in model_run:
        multiplyer=2
    X1 = [[lat_pcb,lon_pcb - .14*multiplyer], [lat_pcb, lon_pcb+.14*multiplyer]]
    X2 = [[lat_pcb-0.08, lon_pcb-0.08], [lat_pcb+0.08, lon_pcb+0.08]]
    X3 = [[lat_pcb+0.08, lon_pcb-0.08], [lat_pcb-0.08, lon_pcb+0.08]]
    
    # take mean of vert motion between lvls 25-48 approx 500m - 1500m
    with warnings.catch_warnings():
        # ignore warning from collapsing non-contiguous dimension:
        warnings.simplefilter('ignore')
        wmean = w[:,25:48,:,:].collapsed('model_level_number', iris.analysis.MEAN)
    h0,h1 = wmean.coord('level_height').bounds[0]
    
    # pull data from masked arrays from cubes
    zi = z.data.data
    topogi = topog.data.data
    # for each timestep:
    for i in range(len(ffdtimes)):
        qci = qc[i].data.data
        wi = w[i].data.data
        ui = u[i].data.data
        wmeani = wmean[i].data.data 
        ffi = None
        if ff is not None:
            ffi = ff[i].data.data
        
        ## First make the left to right figure
        left_right_slice(qci, ui, wi, zi, topogi, lat, lon, X1)
        fio.save_fig(model_run, _sn_, ffdtimes[i], plt, subdir='LR1')
        
        ## second make the full pyrocb plot:
        
        # datetime timestamp for file,title
        stitle = ffdtimes[i].strftime("Vertical motion %Y %b %d %H:%M (UTC)")
        wmeantitle='Mean(%3.0fm - %4.0fm)'%(h0,h1)
        
        pyrocb(w=wi, u=ui, qc=qci, z=zi, wmean=wmeani, topog=topogi,
               lat=lat, lon=lon, transect1=X1, transect2=X2, transect3=X3,
               ff=ffi,
               wmeantitle=wmeantitle,
               extentname=extentname)
        # Save figure into animation folder with numeric identifier
        plt.suptitle(stitle)
        fio.save_fig(model_run, _sn_, ffdtimes[i], plt, dpi=200)

if __name__ == '__main__':
    
    ## check to see where pcb are occurring
    if True:
        sample_showing_grid(model_run="sirivan_run2_hr", extentname='sirivans', HSkip=8)
    

    ## New zoomed, moving pyrocb plotting
    if False:
        waroona_second_half = [datetime(2016,1,5,15)+ timedelta(hours=12+x) for x in range(12)]
        sirivan_check = __PCB_occurrences__['sirivan_run1']['time']
        #moving_pyrocb(model_run='waroona_run3', hours=waroona_second_half)
        moving_pyrocb(model_run='sirivan_run1', 
                      #hours=sirivan_check,
                      xlen=0.3)
    
    ## Run sample for waroona_run2
    if False:
        for hour in waroona_second_half:
            pyrocb_model_run('waroona_run3', dtime=hour)
    
    ### These are the first pyrocb plots I made (3 transects, not moving)
    if False:
        model_runs = ['waroona_run3','waroona_run2','sirivan_run1']
        for mr in model_runs :
            dtimes = fio.model_outputs[mr]['filedates']
            if testing:
                dtimes = dtimes[0:2]
            for dtime in dtimes:
                pyrocb_model_run(model_run=mr, dtime=dtime)
