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

# transects showing the pyroCB in model extents
# waroona (old run) gets pyrocb around Jan 6, 0600-0830
# Sirivan gets some pyrocb around Feb 12, 0600 - 0800 UTC
__PYRO_TRANSECTS__ = {'waroona':{'X1':[[-32.89, 115.90], [-32.88, 116.19]],
                                 'X2':[[-32.95, 116.09], [-32.78, 116.16]],
                                 'X3':[[-32.95, 116.15], [-32.78, 116.09]],
                                 'LR1':[[-32.87, 115.90], [-32.87, 116.18]]}, # left to right transect, to allow easy quiver on (u,w) dims
                      'sirivan':{'X1':[[-32.15, 149.40], [-32.00, 150.10]],
                                 'X2':[[-32.20, 149.60], [-31.85, 149.85]],
                                 'X3':[[-32.15, 149.90], [-31.85, 149.60]],
                                 'LR1':[[-32.50, 149.50], [-32.50, 150.00]]}}

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

def sample_showing_grid(model_run='waroona_run3'):
    """
    Show each hour the latlon grid and vertical motion at 4000,5000,6000 metres
    """
    all_hours = fio.model_outputs[model_run]['filedates']
    extentname=model_run.split('_')[0]
    extent=plotting._extents_[extentname]
    
    for di,hour in enumerate(all_hours):
        ## Read hour of output
        try:
            cubes=fio.read_model_run(model_run,fdtime=hour,extent=extent,
                                     add_winds=True)
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
        FF, = fio.read_fire(model_run, dtimes=[hour], extent=extent)
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
        

def pyrocb(w, qc, z, wmean, topog, lat, lon,
           transect1, transect2, transect3,
           ff=None,
           wmeantitle='Average vertical motion',
           extentname=None,
           ztop=15000,
           cloud_threshold=constants.cloud_threshold,
           ext='.png',
           dpi=300,
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
    wslice, xslice, zslice = plotting.transect_w(w, z, lat, lon, start, end,
                                                 npoints=100, title='',
                                                 topog=topog, ztop=ztop,
                                                 lines=None, colorbar=False)
    ## add cloud outlines
    ## Add contour where clouds occur
    qcslice = utils.cross_section(qc, lat, lon, start, end, npoints=100)
    with warnings.catch_warnings():
        # ignore warning when there are no clouds:
        warnings.simplefilter('ignore')
        plt.contour(xslice,zslice,qcslice,np.array([cloud_threshold]),colors='k')

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

def moving_pyrocb(model_run='waroona_run2',subset=True):
    """
    follow pyrocb with a transect showing vert motion and pot temp
    set subset to false if running on full dataset (on NCI)
    """
    extentname=model_run.split('_')[0]
    extent = plotting._extents_[extentname]
    clevs_vertwind = np.union1d(np.union1d(2.0**np.arange(-2,6),
                                           -1*(2.0**np.arange(-2,6))),
                                np.array([0]))
    ztop=14000
    cloud_threshold = constants.cloud_threshold
    
    def lin_space_transect(start,end,steps):
        j0,k0,l0,m0 = start
        j1,k1,l1,m1 = end
        out = [ np.linspace(a,b,steps) for [a,b] in [[j0,j1],[k0,k1],[l0,l1],[m0,m1]] ]
        outt= np.array(out).T # transpose
        loutt = list([list(intt) for intt in outt])
        return loutt
        
    # zoomed in transects for pyrocb:
    tran = {'sirivan_run1':[[-32.05,149.5,-32.11,149.8],]*9+
                           [[-32.07,149.6,-32.15,150.0],]*5 + # up to 0331utc
                           [[-32.05, 149.6, -32.1, 150.0], #0401
                            [-32.038, 149.58, -32.075, 150.0], # 0431
                            [-32.026, 149.56, -32.05, 150.0], # 0501
                            [-32.025, 149.55, -32.04, 150.0], # 0531
                            [-32.025, 149.65, -32.04, 150.1], # 0601
                            [-32.00, 149.53, -32.03, 150.0], # 0631
                            [-31.98, 149.58, -32.02, 150.1],  # 0701
                            [-31.96, 149.53, -31.92, 150.1],  # 0730
                            [-31.96, 149.55, -31.88, 150.1]]+  # 0801
                           lin_space_transect([-31.94,149.5,-31.90,150.1],
                                              [-31.83,149.5,-31.82,150.2],25),
            'waroona_old':[[-32.9, 116.15, -32.86, 116.19],]*15 + # up to 2201
                           [[-32.9, 116.15, -32.86, 116.19], # 2231
                            [-32.898, 116.144, -32.86, 116.19], # 2301
                            [-32.896, 116.138, -32.86, 116.19], # 2330
                            [-32.893, 116.133, -32.86, 116.19], # 0000
                            [-32.891, 116.127, -32.86, 116.19], # 0030
                            [-32.888, 116.122, -32.86, 116.19], # 0100
                            [-32.886, 116.116, -32.86, 116.19], # 0130
                            [-32.884, 116.111, -32.86, 116.19], # 0200
                            [-32.882, 116.105, -32.86, 116.19],]+ # 0230
                            lin_space_transect([-32.88, 116.1, -32.86, 116.19],
                                               [-32.92, 116.1, -32.84, 116.17],13), # 0300 - 0830
            # Run 1 starts later (8 hrs later!?) starts at 0030 UTC
            'waroona_run1':[[-32.89, 116.15, -32.886, 116.19],]*60 + \
                           [[-32.888, 116.14, -32.886, 116.19],]*3 + \
                           [[-32.886, 116.13, -32.886, 116.19],]*3 + \
                           [[-32.886, 116.12, -32.886, 116.19],]*3 + \
                           [[-32.886, 116.11, -32.886, 116.19],]*3 + \
                           lin_space_transect([-32.886, 116.11, -32.886, 116.19],
                                              [-32.881, 115.95, -32.885, 116.09],57) + \
                           lin_space_transect([-32.879, 115.95, -32.885, 116.09],
                                              [-32.872, 115.92, -32.885, 116.09],15), #1240 - 1440
            'waroona_run2':[[-32.89, 116.10, -32.886, 116.19],]*60 + \
                           lin_space_transect([-32.89, 116.09, -32.886, 116.19],
                                              [-32.872, 115.85, -32.885, 116.07],84),
           }
    if subset:
        # subset has 2 outputs/hr, full dataset has 6
        tran['waroona_run1'] = tran['waroona_run1'][::3]
        tran['waroona_run2'] = tran['waroona_run2'][::3]
    tran['waroona_run2uc'] = tran['waroona_run2']
        
    transects = tran[model_run]
    #datetimes = fio.model_outputs[model_run]['filedates']
    ## read um output over extent [t, lev, lat, lon]
    cubes = fio.read_model_run(model_run, extent=extent, add_topog=True)
                               #add_z=True, add_winds=True, add_topog=True)
    #print("DEBUG:",cubes)
    w, = cubes.extract('upward_air_velocity')
    #u, = cubes.extract('u')
    qc, = cubes.extract('qc')
    topog0, = cubes.extract('surface_altitude')
    topog = topog0.data
    #z, = cubes.extract('z_th') 
    
    lat = w.coord('latitude').points
    lon = w.coord('longitude').points

    ## fire front
    ffdtimes = utils.dates_from_iris(w)
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
        wi = w[i].data.data
        ## calc zth
        # repeat surface pressure along z axis
        reppmsl = np.repeat(pmsl[i].data[np.newaxis,:,:],nz, axis=0)
        zth = -(287*300/9.8)*np.log(p[i].data/reppmsl)
        
        a,b,c,d = transect
        start=[a,b]
        end=[c,d]
        
        ## map showing transect
        fig = plt.figure(figsize=[10,16])
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
        
        wslice, xslice, zslice = plotting.transect_w(wi, zth,
                                                     lat, lon, start, end,
                                                     npoints=100, title='',
                                                     topog=topog, ztop=ztop,
                                                     contours=clevs_vertwind,
                                                     lines=None, colorbar=True)
        ## add cloud outlines
        ## Add contour where clouds occur
        qcslice = utils.cross_section(qci, lat, lon, start, end, npoints=100)
        with warnings.catch_warnings():
            # ignore warning when there are no clouds:
            warnings.simplefilter('ignore')
            plt.contour(xslice,zslice,qcslice,np.array([cloud_threshold]),colors='k')
    
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
    X1 = __PYRO_TRANSECTS__[extentname]['X1']
    X2 = __PYRO_TRANSECTS__[extentname]['X2']
    X3 = __PYRO_TRANSECTS__[extentname]['X3']
    LR1 = __PYRO_TRANSECTS__[extentname]['LR1']
    
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
        left_right_slice(qci, ui, wi, zi, topogi, lat, lon, LR1)
        fio.save_fig(model_run, _sn_, ffdtimes[i], plt, subdir='LR1')
        
        ## second make the full pyrocb plot:
        
        # datetime timestamp for file,title
        stitle = ffdtimes[i].strftime("Vertical motion %Y %b %d %H:%M (UTC)")
        wmeantitle='Mean(%3.0fm - %4.0fm)'%(h0,h1)
        
        pyrocb(w=wi, qc=qci, z=zi, wmean=wmeani, topog=topogi,
               lat=lat, lon=lon, transect1=X1, transect2=X2, transect3=X3,
               ff=ffi,
               wmeantitle=wmeantitle,
               extentname=extentname)
        # Save figure into animation folder with numeric identifier
        plt.suptitle(stitle)
        fio.save_fig(model_run, _sn_, ffdtimes[i], plt, dpi=200)



if __name__ == '__main__':
    
    model_runs = ['waroona_run1','waroona_old','sirivan_run1']
    testing=False
    
    ## New zoomed, moving pyrocb plotting
    moving_pyrocb(model_run='waroona_run2uc',subset=True)
    
    ## Run sample for waroona_run2
    #pyrocb_model_run('waroona_run2', dtime=datetime(2016,1,5,15))
    
    ### These are the first pyrocb plots I made (3 transects, not moving)
    #for mr in model_runs :
    #    dtimes = fio.model_outputs[mr]['filedates']
    #    if testing:
    #        dtimes = dtimes[0:2]
    #    for dtime in dtimes:
    #        pyrocb_model_run(model_run=mr, dtime=dtime)
