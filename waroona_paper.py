#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 16:28:18 2019
    Show fire spread and intensity over time
@author: jesse
"""


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects, colors, image, ticker, patches
import matplotlib.dates as mdates
from matplotlib.ticker import FormatStrFormatter,LinearLocator
import matplotlib.patheffects as PathEffects
import numpy as np
# for legend creation:
from matplotlib.lines import Line2D
from datetime import datetime,timedelta

import iris
import string
import warnings

from cartopy import crs as ccrs

from utilities import plotting, utils, fio, constants
from pyrocb import map_with_transect, pcb_occurrences, transect_plus_stream
from emberstorm import _emberstorm_centres_, emberstorm_centres, add_vertical_contours

###
## GLOBALS
###
_sn_ = 'waroona_paper'

## plotting defaults
matplotlib.rcParams['font.size'] = 16.0
matplotlib.rcParams["text.usetex"]  = False     # I forget what this is for, maybe allows latex labels?
matplotlib.rcParams["legend.numpoints"] = 1         # one point for marker legends
matplotlib.rcParams["figure.figsize"]   = (10, 9)    # Default figure size (xinches,yinches)
matplotlib.rcParams["axes.titlesize"]   = 18        # title font size
matplotlib.rcParams["figure.titlesize"] = 20        # figure suptitle size
matplotlib.rcParams["axes.labelsize"]   = 16        #
matplotlib.rcParams["xtick.labelsize"]  = 14        #
matplotlib.rcParams["ytick.labelsize"]  = 14        #
matplotlib.rcParams['image.cmap'] = 'plasma'        # Colormap default
matplotlib.rcParams['axes.formatter.useoffset'] = False # another one I've forgotten the purpose of
matplotlib.rcParams['figure.dpi'] = 600 # important for literature quality

def add_labels(fig,skips=None):
    """
    """
    ax_list = fig.axes
    #print("DEBUG: add_labels:",np.shape(ax_list),ax_list)
    letters=string.ascii_uppercase[:len(ax_list)]
    for i, label in enumerate(letters):
        j=0
        if skips is not None:
            if i in skips:
                j = j+1
        if i+j >= len(ax_list):
            break
        ax=ax_list[i+j]
        ax.text(0.075, 1.075, label, transform=ax.transAxes,
           fontsize=16, 
           #fontweight='bold',
           va='top', 
           ha='right')

def fireplan(ff, fire_contour_map = 'autumn',
             show_cbar=True, cbar_XYWH= [0.65, 0.63, .2, .02],
             extentname=None,
             color_by_day=None,
             last_hour=None,
             **kwtiffargs):
    '''
    show satellite map of extent, with fire front over time overplotted

    ARGUMENTS:
        ff: iris.cube.Cube with time, longitude, latitude dimensions
        fire_contour_map: how will firefront contours be coloured
        show_cbar: bool
            draw a little colour bar showing date range of contours?
        cbar_XYHW: 4-length list where to put cbar
        color_by_day: dictionary for colouring by day
            {"6":"red",...,"8":None}
            Set contour colours based on "%-d" of localtime
        fig,... arguments to plotting.map_tiff_qgis()
    '''
    lon,lat = ff.coord('longitude').points, ff.coord('latitude').points
    _,ny,nx = ff.shape
    extent = [np.min(lon),np.max(lon), np.min(lat),np.max(lat)]

    # Get datetimes from firefront cube
    ftimes = utils.dates_from_iris(ff)
    
    # How many fire contours do we have?
    minff = np.min(ff.data,axis=(1,2))
    subset_with_fire = minff < 0
    ff_f = ff[subset_with_fire]
    ftimes = ftimes[subset_with_fire]
    
    # just read hourly by checking when the hour changes
    hourinds = [ftimes[i].hour != ftimes[i-1].hour for i in range(len(ftimes)+1)[:-1]]
    #hourinds = [(ft.minute==0) and (ft.second==0) for ft in ftimes]
    nt = np.sum(hourinds)
    
    ## fire contour colour map
    cmap = matplotlib.cm.get_cmap(fire_contour_map)
    rgba = cmap(np.linspace(0,1,nt))

    ## PLOTTING

    # First show satellite image and locations
    extentname='waroona'
    for k,v in kwtiffargs.items():
        print("    :",k,v)
    fig, ax = plotting.map_tiff_qgis(
        fname=extentname+'.tiff',
        extent=extent, 
        **kwtiffargs,
        )
    plotting.map_add_locations_extent(extentname, hide_text=False, nice=True)
    plt.xlabel('latitude')
    plt.ylabel('longitude')
    # plot contours at each hour
    tstamp=[]
    for ii,dt in enumerate(ftimes[hourinds]):
        if last_hour is not None:
            if dt > last_hour+timedelta(minutes=5):
                break
        if 'sirivan' in extentname:
            mr='sirivan_run3'
        elif 'waroona' in extentname:
            mr='waroona_run3'
        else:
            mr='KI_run0'
        offset=fio.run_info[mr]['UTC_offset']
        LT = dt + timedelta(hours=offset)
        color=[rgba[ii]]
        # Colour waroona by day of firefront
        if color_by_day is not None:
            dnum = LT.strftime("%-d")
            color = color_by_day[dnum]
            if color is None:
                continue # Skip if None used as colour
        tstamp.append(LT.strftime('%b %d, %H%m(local)'))

        ffdata = ff_f[hourinds][ii].data.data
        
        linewidth=1 + (ii == nt-1) # make final hour thicker
        plotting.map_fire(ffdata,lat,lon,
                          linewidths=linewidth,
                          colors=color)
        #print("INFO:Contour for ",LT,"(local time) plotted")
    if (color is not None) and (last_hour is None):
        # final contour gets bonus blue line
        final_line=plt.contour(lon,lat,ff_f[-1].data.data, np.array([0]), 
                               linestyles='dotted',
                               colors='cyan', linewidths=1)
        clbls = plt.clabel(final_line,[0],fmt=LT.strftime('%H:%M'), 
                           inline=True, colors='wheat')
        plt.setp(clbls, path_effects=[patheffects.withStroke(linewidth=3, foreground="k")])

    ## Add tiny colour bar showing overall time of fire
    if show_cbar:
        # create an axis somewhere to add a colorbar
        cax = fig.add_axes(cbar_XYWH)
        #cax.axes.get_xaxis().set_visible(False)
        cax.axes.get_yaxis().set_visible(False) # horizontal alignment, no yax
        cbar = matplotlib.colorbar.ColorbarBase(cax,cmap=plt.get_cmap(fire_contour_map), orientation='horizontal')
        cbar.set_ticks([0,1])
        cbar.set_ticklabels([tstamp[0],tstamp[-1]])
        # make xtick labels have black outline on wheat text color
        cbxtick_obj = plt.getp(cax, 'xticklabels') # get xtick labels
        plt.setp(cbxtick_obj, color='wheat',
                 path_effects=[patheffects.withStroke(linewidth=3, foreground="k")])
        # return focus to newly created plot
        plt.sca(ax)

    return fig, ax


def plot_X_transect(w, u, qc, z, wmean, topog, lat, lon,
           transect1, transect2, transect3,
           ff=None, sh=None,
           wmeantitle='Average vertical motion',
           extentname=None,
           show_windstream=False,
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
    #plotting.init_plots()
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
                              cbargs={}, # effectively turns off cbar
                              levels = clevs_vertwind, 
                              cmap=plotting._cmaps_['verticalvelocity'],
                              norm=colors.SymLogNorm(0.25,base=2.), 
                              #cbarform=ticker.ScalarFormatter(), 
                              #clabel='m/s',
                              )
    plt.title(wmeantitle)
    
    # add nearby towns
    if extentname is not None:
        plotting.map_add_locations_extent(extentname,nice=True, fontsize=17)
    
    ### transect plots
    ###
    ax2=plt.subplot(3,1,2)
    
    transect_plus_stream(w,u,qc,topog,z,lat,lon,start,end,ztop,
                         contours=clevs_vertwind,
                         showcolorbar=False, # making our own colorbar
                         show_windstream=show_windstream,
                         ) 
    
    # REVIEW: changed to km forccing labels
    plt.ylabel('Alt (km)')
    plt.yticks([0.,5000.,10000.,15000.],[0,5,10,15])
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
                                                      sh=sh,
                                                      topog=topog, 
                                                      colorbar=False,
                                                      ztop=ztop,
                                                      lines=None)
        ## add cloud outlines
        ## Add contour where clouds occur
        qcslicex,_,_ = utils.transect(qc, lat, lon, startx, endx, nx=100)
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

    plt.sca(ax31)
    plt.yticks([0.,5000.,10000.,15000.],[0,5,10,15])
        
    # add colour bar
    axes=[0.85, 0.20, 0.03, 0.6] 
    #f.subplots_adjust(top=0.95)
    #f.subplots_adjust(hspace=0.05)
    f.subplots_adjust(right=axes[0]-0.01)
    cbar_ax = f.add_axes(axes)
    cb = f.colorbar(cs, cax=cbar_ax, 
                    format=ticker.ScalarFormatter())
    # -ve labelpad moves label closer to colourbar
    cb.set_label('ms$^{-1}$', labelpad=-4)


def topdown_emberstorm(fig=None, subplot_row_col_n=None, ax=None,
                       extent=[115.68, 116.15, -33.025,-32.79],
                       lats=None,lons=None, ff=None, sh=None, 
                       u10=None, v10=None, 
                       wmap=None, wmap_height=None,
                       topog=None,
                       annotate=True, showlatlons=True,
                       sh_kwargs={},
                       ):
    """
    Top down view of Waroona/Yarloop, adding fire front and heat flux and 10m winds
    ARGUMENTS:
        ax: plotting axis, if this is provided then no backdrop is drawn (assume axis already has backdrop)
            In this case just winds/fire/etc will be overplotted
        lats/lons are 1D arrays, required if adding other stuff
        ff [lats,lons] : firefront array
        sh [lats,lons] : sensible heat flux
        u10 [lats,lons] : 10m altitude longitudinal winds
        v10 [lats,lons] : 10m altitude latitudinal winds
        topog [lats,lons] : surface altitude - can use this instead of tiff
        topog_contours : list of values to contour topography at - default 50m
        annotate: True if annotations are desired for winds and heat flux
        showlatlons: True if latlons should be added to border
    RETURNS:
        fig, ax
    """
    annotate_font_size=13
    # if we already have an axis, assume the backdrop is provided
    if ax is None:
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
            
            ## Add waroona, yarloop if possible
            for txt in ['Waroona','Yarloop']:
                ax.annotate(txt, xy=np.array(constants.latlons[str.lower(txt)])[::-1],
                            xycoords="data", # lat lon xy as coords are platecarree
                            fontsize=14, ha="center",
                            color='k',
                            path_effects=[PathEffects.withStroke(linewidth=2,
                                                                 foreground="w")])
            
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    
    if ff is not None:
        # add firefront
        cs_ff = plotting.map_fire(ff,lats,lons,colors=['orange'],linewidths=[2])
    if sh is not None:
        # add hot spots for heat flux
        # default kwargs for sh plot
        if 'alpha' not in sh_kwargs:
            sh_kwargs['alpha']=1.0
        if 'cbar_kwargs' not in sh_kwargs:
            sh_kwargs['cbar_kwargs'] = {'label':"Wm$^{-2}$"}
        cs_sh, cb_sh = plotting.map_sensibleheat(sh,lats,lons,colorbar=False,**sh_kwargs)
        if annotate:
            plt.annotate(text="max heat flux = %6.1e W/m2"%np.max(sh),
                         xy=[0,1.06],
                         xycoords='axes fraction', 
                         fontsize=annotate_font_size)
    if u10 is not None:
        # winds, assume v10 is also not None
        s10 = np.hypot(u10,v10)
        speedmax=20 # what speed for thickest wind streams
        lwmax_winds=5 # how thick can the wind streams become
        lw10 = utils.wind_speed_to_linewidth(s10, lwmax=lwmax_winds, speedmax=speedmax)
        # higher density if using topography instead of OSM
        density=(0.6,0.5) if topog is None else (0.75,0.7)
        # REVIEW: UPDATE removed streamplot
        #plt.streamplot(lons,lats,u10,v10, 
        #               linewidth=lw10, 
        #               color='k',
        #               density=density,
        #               )
        
        if annotate:
            #plt.annotate("10m wind linewidth increases up to %dms$^{-1}$"%(speedmax),
            #             xy=[0,1.09], 
            #             xycoords="axes fraction", 
            #             fontsize=10)
            plotting.annotate_max_winds(s10, text="10m wind max = %5.1f m/s",
                                        xytext=[0,1.025],
                                        fontsize=annotate_font_size)
    
    if wmap is not None:
        add_vertical_contours(wmap,lats,lons,
                              wmap_height=wmap_height,
                              wmap_levels=[1,3],
                              annotate=True,
                              xy0=[0.73,1.07],
                              annotate_font_size=annotate_font_size)
        
    # set limits back to latlon limits
    ax.set_ylim(ylims[0],ylims[1])
    ax.set_xlim(xlims[0],xlims[1])
    # 115.8, 116.1, -32.92,-32.82
    if showlatlons:
        #xticks=np.arange(115.8,116.11,0.05)
        #plt.xticks(xticks,xticks)
        #yticks=np.arange(-32.92,-32.805,0.03)
        #plt.yticks(yticks,yticks)
        ax.xaxis.set_major_locator(LinearLocator(numticks=5))
        ax.yaxis.set_major_locator(LinearLocator(numticks=5))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    
    
    return fig, ax



def transect_emberstorm(u,v,w,z,
                        lats,lons,
                        transect,
                        topog,
                        ztop=700,
                        sh=None,
                        theta=None,
                        theta_contourargs={},
                        wind_contourargs={},
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
    retdict = {} # return info for whatever use
    
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
    wind_contourargs['ztop']=ztop
    wind_contourargs['npoints']=npoints
    wind_contourargs['topog']=topog
    wind_contourargs['sh'] = sh
    if 'title' not in wind_contourargs:
        wind_contourargs['title']=""
    if 'lines' not in wind_contourargs:
        wind_contourargs['lines']=None
    
    slices, slicex, slicez = plotting.transect_s(s,z, 
                                            lats,lons, 
                                            start, end,
                                            **wind_contourargs)
    
    # save the max windspeed and location
    mlocs = utils.find_max_index_2d(slices)
    retdict['s'] = slices
    retdict['max_s'] = slices[mlocs]
    retdict['max_s_index'] = mlocs
    retdict['x'] = slicex
    retdict['y'] = slicez
    # east-west and vertical winds on transect
    sliceu,_,_ = utils.transect(u,lats,lons,start,end,nx=npoints, z_th=z)
    slicew,_,_ = utils.transect(w,lats,lons,start,end,nx=npoints, z_th=z)
    
    # Streamplot
    # REVIEW: changed colour to white
    plotting.streamplot_regridded(slicex,slicez,sliceu,slicew,
                                  density=(1,1), 
                                  color='white',
                                  zorder=1,
                                  #linewidth=np.hypot(sliceu,slicew), # too hard to see what's going on
                                  minlength=0.8, # longer minimum stream length
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
            theta_contourargs['alpha'] = 1.0
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
        plt.clabel(contours, inline=True, fontsize=13)# REVIEW: upped from 10->13
    return retdict

def weather_series(
        day2=False,
        extent=None,
        mapname=None,
        showfirelinehour=None
        ):
    """
        time series of surface weather
    """
    model_run='waroona_run3'
    test=False
    HSkip=None
    showFP=False
    showPFT=False
    showwsmax=False
    showmap=True
    showfirelines=(not day2)
    showQC=True
    localtime=True
                   
    
    extentname=model_run.split('_')[0]
    # zoomed extent for analysis
    extentnamez = extentname + 'z'
    if extent is None:
        extent=constants.extents[extentnamez]
    
    qc_thresh = constants.cloud_threshold
    
    # READ EVERYTHING, SUBSET, CALC DESIRED METRICS, PLOT METRICS
    modeltimes = fio.run_info[model_run]['filedates']
    fdtimes = modeltimes[-24:] if day2 else modeltimes[:24]

    cubes = fio.read_model_run(model_run, extent=extent,
                               HSkip=HSkip,
                               fdtime=fdtimes)
    ctimes = utils.dates_from_iris(cubes[0])
    ## read all the fire data, and lats/lons
    ## 10m u and v winds (already destaggered) are in fire model output
    ff,sh,u10, v10 = fio.read_fire(model_run=model_run, dtimes=ctimes,
                           extent=extent, HSkip=HSkip,
                           firefront=True, sensibleheat=True,
                           day2=day2, day1=(not day2),
                           wind=True)
    
    lons,lats = sh.coord('longitude').points, sh.coord('latitude').points
    assert sh is not None, "Missing sensible heat file"
    # datetimes list, convert to local time if desired
    ftimes = utils.dates_from_iris(sh)
    ftimes_lt = ftimes
    ctimes_lt = ctimes
    offset_hrs = fio.run_info[model_run]['UTC_offset']
    if localtime:
        offset = timedelta(hours=offset_hrs)
        ftimes_lt = np.array([ft + offset for ft in ftimes ])
        ctimes_lt = np.array([ct + offset for ct in ctimes ])
    
    # spatial sum over time of sensible heat gives firepower
    # gives GWatts
    firepower = np.sum(utils.firepower_from_cube(sh), axis=(1,2))
    
    # remove zeros:
    prefire = np.isclose(np.cumsum(firepower), 0)
    firepower[prefire] = np.NaN
    
            
    ## Read PFT
    # calculated using kevin's code, stored as GWatts
    lat,lon = plotting._latlons_["fire_%s_upwind"%extentname]
    if showPFT:
        pft, ptimes, _, _ = fio.read_pft(model_run,lats=lat,lons=lon)
        ptimes_lt = ptimes
        #pft = np.nanmean(pft,axis=(1,2)) # average spatially
        if localtime:
            ptimes_lt = np.array([pt + offset for pt in ptimes ])
    
    ## get temperature, RH, cloud
    qTqc = cubes.extract(['specific_humidity','air_temperature','qc'])
    q,T,qc = qTqc

    # closest to 500m height index
    index_500m = np.argmin(np.abs(utils.height_from_iris(q)-500))
    assert index_500m > 0, "Index500 didn't work = %d"%index_500m
    clats, clons = q.coord('latitude').points, q.coord('longitude').points
    # compute RH from specific and T in kelvin
    T.convert_units('K')
    # just want surface for Temp and RH
    q = q[:,0,:,:]
    T = T[:,0,:,:].data
    # qc stuff takes heaps of RAM
    #qc_weight[np.isnan(qc_weight)] = 0.0 # return the zeros
    RH = utils.relative_humidity_from_specific(q.data, T)
    RH = np.mean(RH, axis=(1,2))
    T = np.mean(T, axis=(1,2))
    ## get wind speed/ wind dir
    u, v = cubes.extract(['x_wind','y_wind'])
    ## wind speed at 10m is output by the fire model
    u10=u10.data #np.swapaxes(u10.data, 1,2) # time, lon, lat -> time, lat, lon
    v10=v10.data #np.swapaxes(v10.data, 1,2) # also for v10
    u500=u[:,index_500m,:,:]
    v500=v[:,index_500m,:,:]
    
    # DESTAGGER u and v using iris interpolate
    u500 = u500.interpolate([('longitude',v.coord('longitude').points)],
                        iris.analysis.Linear())
    v500 = v500.interpolate([('latitude',u.coord('latitude').points)],
                        iris.analysis.Linear())
    ws10_all = utils.wind_speed_from_uv(u10,v10)
    ws10_q1 = np.quantile(ws10_all, 0.25, axis=(1,2))
    ws10_q3 = np.quantile(ws10_all, 0.75, axis=(1,2))
    ws10_max = np.max(ws10_all, axis=(1,2))
    ws10 = np.mean(ws10_all, axis=(1,2))
    ws500_all = utils.wind_speed_from_uv_cubes(u500,v500).data
    ws500_q1 = np.quantile(ws500_all, 0.25, axis=(1,2))
    ws500_q3 = np.quantile(ws500_all, 0.75, axis=(1,2))
    #ws500_max = np.max(ws500_all, axis=(1,2))
    ws500 = np.mean(ws500_all, axis=(1,2))
    # mean wind direction based on mean u,v
    u10_mean = np.mean(u10,axis=(1,2))
    v10_mean = np.mean(v10,axis=(1,2))
    u500_mean = np.mean(u500.data,axis=(1,2)) 
    v500_mean = np.mean(v500.data,axis=(1,2))
    #wd0 = utils.wind_dir_from_uv(u10_mean,v10_mean)
    
    #######################
    #### PLOTTING PART ####
    #######################
    
    # figure and multiple y axes for time series
    fig = plt.figure(figsize=[10,11])
    dx,dy = extent[1]-extent[0], extent[3]-extent[2]
    extentplus = np.array(extent)+np.array([-0.3*dx,0.3*dx,-0.2*dy,0.2*dy])
    if showmap:
        # map with extent shown
        if mapname is None:
            mapname = "%s.tiff"%extentname
        fig, ax_map = plotting.map_tiff_qgis(fname=mapname, 
                                             fig=fig,
                                             extent=list(extentplus), 
                                             subplot_row_col_n = [3,1,1],
                                             #locnames=[extentname,]
                                             )
        # REVIEW: may try this
        if day2:
            ax_map.set_aspect('equal')
            plotting.map_add_locations_extent('waroona_day2',
                                              hide_text=False,
                                              nice=True,
                                              abbr=True, #REVIEW: added changing name to just letter
                                              fontsize=16,
                                              )
        else:
            plotting.map_add_locations_extent(extentnamez,
                                              hide_text=False,
                                              nice=True,
                                              abbr=True, #REVIEW:...
                                              fontsize=16,
                                              )
        # Add rectangle
        botleft = extent[0],extent[2]
        ax_map.add_patch(patches.Rectangle(xy=botleft,
                                           width=dx,
                                           height=dy,
                                           fill=False,
                                           edgecolor='blue',
                                           linewidth=3, #REVIEW: updated 2->3
                                           alpha=0.6, 
                                           ))
        
        # add fire outline
        plt.sca(ax_map)
        # Firefront shown at each hour increment:
        if showfirelines:
            # look at roughly every 2 hours, want final outline too
            for ffhind in range(0+len(ctimes)//12-1,len(ctimes),len(ctimes)//12):
            #for ffhour in fdtimes[3::4]: 
                ffline = plotting.map_fire(ff[ffhind].data, lats, lons,
                        alpha=0.7, #REVIEW: updated 0.4->0.7
                        linewidths= 2.0, #REVIEW: added at 2.0 
                        colors=['orange'], #REVIEW: added for jeff colorblind
                        )
        if showfirelinehour is not None:
            ffsample,= fio.read_fire(model_run=model_run, dtimes=[showfirelinehour],
                           extent=list(extentplus), HSkip=HSkip,
                           firefront=True, 
                           day2=day2, day1=(not day2),)
            ffsample_lats = ffsample.coord('latitude').points
            ffsample_lons = ffsample.coord('longitude').points
            ffline = plotting.map_fire(ffsample[0].data, ffsample_lats,ffsample_lons,
                    alpha=0.6, #REVIEW: updated 0.4->0.6
                    linewidths=2.0, #REVIEW: added at 2.0,
                    colors=['orange'], #REVIEW: added for jeff colorblind
                    )

        ax_map.set_title("Averaged area")
    
    ax_fp = plt.subplot(2+showmap,1,1+showmap) # firepower axis
    color_fp = "orange"
    ax_T = ax_fp.twinx() # Temperature
    color_T = "red"
    if (not showFP) and (not showPFT):
        ax_RH = ax_fp # just use left axis if no firepower
    else:
        ax_RH = ax_fp.twinx() # Rel Humid
        # offsets for extra y axes, and make just the one spine visible
        ax_RH.spines["right"].set_position(("axes", 1.07))
        plotting.make_patch_spines_invisible(ax_RH)
        ax_RH.spines["right"].set_visible(True)
        ax_fp.set_ylabel("GWatts")
        ax_fp.yaxis.label.set_color(color_fp)
    color_RH = "magenta"
    ax_ws = plt.subplot(2+showmap,1,2+showmap, sharex=ax_fp) # Wind speed
    color_ws = 'black'
    color_qc = "darkblue"
    color_wd = "brown"
    # wind speed at 500m on same axis as at 10m
    #ax_ws500 = ax_ws.twinx() # wind speed at 500m
    color_ws500 = 'chocolate'
    
    
    ## Plot firepower, and PFT on one axis
    plt.sca(ax_fp)
    if showPFT:
        line_pft, = plt.plot_date(ptimes_lt, pft, '--', 
                                color=color_fp, label='PFT', alpha=0.6)
    if showFP:
        line_fp, = plt.plot_date(ftimes_lt, firepower, '-',
                                 color=color_fp, label='firepower')
    else:
        ignition=ftimes_lt[~prefire][0]
        hint="" if day2 else "Fire"
        plt.annotate(hint,(mdates.date2num(ignition), 0), xytext=(-60, -5), 
                     textcoords='offset points', arrowprops=dict(arrowstyle='-|>'), 
                     color='red')
    ## plot temperature on right axis
    plt.sca(ax_T)
    line_T, = plt.plot_date(ctimes_lt, T-273.15, '-',
                           color=color_T, label="Temperature",)
    ## plot RH on further right axis
    plt.sca(ax_RH)
    line_RH, = plt.plot_date(ctimes_lt, RH, '-',
                             color=color_RH, label="RH",)

    if showQC:
        qc_sum = qc.collapsed('model_level_number', iris.analysis.SUM)
        qc_frac = np.sum(qc_sum.data > qc_thresh, axis=(1,2))/(len(clats)*len(clons))
        qc_weight = qc_sum.collapsed(['longitude','latitude'], iris.analysis.SUM).data
        #qc_weight[qc_weight <= 0.0001] = np.NaN # take out zeros
        qc_q3 = np.nanpercentile(qc_weight, 75)
        qc_q2 = np.nanpercentile(qc_weight, 50)
        qc_heavy = qc_weight > qc_q3
        qc_mid = (qc_weight > qc_q2) * (qc_weight < qc_q3)
        plt.plot_date(ctimes_lt, qc_frac, 'o',
                      color=color_qc, 
                      fillstyle='none', 
                      mec=color_qc,
                      mew=1,)
        line_qc, = plt.plot_date(ctimes_lt[qc_heavy], qc_frac[qc_heavy], 'o',
                                 color=color_qc,
                                 label="Clouds",
                                 fillstyle='full')
        plt.plot_date(ctimes_lt[qc_mid], qc_frac[qc_mid], 'o',
                      color=color_qc,
                      fillstyle='bottom')
    
    ## Wind speed, wind speed stdev, wind dir quiver
    plt.sca(ax_ws)
    line_ws, = plt.plot_date(ctimes_lt, ws10, 'o',
                             color=color_ws, label="10m wind speed (+IQR)",)
    if showwsmax:
        line_wsmax, = plt.plot_date(ctimes_lt, ws10_max, '^', 
                                    color=color_ws, label="10m max wind speed")
    
    # add quartiles
    plt.fill_between(ctimes_lt, ws10_q3, ws10_q1, alpha=0.6, color='grey')
        
    ## wind direction quiver
    # normalize windspeed for unit length quivers
    wdnorm = np.sqrt(u10_mean**2 + v10_mean**2)
    # dont show quiver at every single point
    qskip = len(wdnorm)//24 # just show 24 arrows
    plt.quiver(ctimes_lt[::qskip], ws10_q1[::qskip], 
               u10_mean[::qskip]/wdnorm[::qskip], 
               v10_mean[::qskip]/wdnorm[::qskip], 
               pivot='mid', 
               alpha=.9, 
               color=color_ws,
               headwidth=3, headlength=2, headaxislength=2,
               width=.004)
    ## Same again but for 500m wind
    #plt.sca(ax_ws500)
    line_ws500, = plt.plot_date(ctimes_lt, ws500, '.',
                                color=color_ws500, 
                                label="500m wind speed (+IQR)",)
    # add quartiles
    plt.fill_between(ctimes_lt, ws500_q3, ws500_q1, alpha=0.4, color=color_ws500)
    
    # normalize windspeed for unit length quivers
    wdnorm500 = np.sqrt(u500_mean**2 + v500_mean**2)
    if not test:
            
        plt.quiver(ctimes_lt[::qskip], ws500_q3[::qskip], 
                   u500_mean[::qskip]/wdnorm500[::qskip], 
                   v500_mean[::qskip]/wdnorm500[::qskip], 
                   pivot='mid', 
                   alpha=.7, 
                   color=color_wd,
                   headwidth=3, headlength=2, headaxislength=2,
                   width=.003)
    # top row:
    
    ax_T.set_ylabel("Temperature (C)")
    ax_T.yaxis.label.set_color(color_T)
    ax_RH.set_ylabel("RH and cloud (frac)")
    ax_RH.yaxis.label.set_color(color_RH)
    
    # bottom row:
    ax_ws.set_ylabel("wind speed (m/s)")
    ax_ws.yaxis.label.set_color(color_ws)
    ax_ws.set_xlabel("UTC + %d"%offset_hrs)
    #ax_ws500.set_ylabel("500m wind speed (m/s)")
    #ax_ws500.yaxis.label.set_color(color_ws500)
    
    ## Plot periphery
    lines = [line_T, line_RH,]
    if showQC:
        lines.append(line_qc)
    if showFP:
        lines.append(line_fp)
    if showPFT:
        lines.append(line_pft)
    labels = [l.get_label() for l in lines]
    #REVIEW: moved legend for day2
    if day2:
        legloc='center right'
    else:
        legloc='center left'
    ax_T.legend(lines, labels, loc=legloc)
    #plt.suptitle(model_run,y=0.9)
    
    if showwsmax:
        lines2 = [line_ws, line_wsmax, line_ws500]
    else:
        lines2 = [line_ws, line_ws500]
    labels2 = [l.get_label() for l in lines2]
    ax_ws.legend(lines2,labels2, loc='best')
    
    # format x-ticks date times
    ax_fp.xaxis.set_major_locator(mdates.HourLocator(interval=3))
    ax_fp.xaxis.set_minor_locator(mdates.HourLocator())
    ax_fp.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    
    fig.tight_layout(rect=[0,0,.99,1]) # left, bottom, right, top

    pname = 'fig3.png' if day2 else 'fig2.png'
    #add_labels(fig,skips=[1,]) # doesn't work well on figures with twinned axis
    fio.save_fig(model_run, _sn_, pname, plt=plt,dpi=600)


def vert_motion_slices(qc,w,lh,lat,lon,
                       ff=None,
                       extentname=None,
                       cloud_threshold=constants.cloud_threshold,
                       three_by_three=False,
                       ):
    '''
    44i showing vertical motion contourf plots, at different model levels
    Trying to see pyroCB cloud
    ARGUMENTS:
        qc [lev,lat,lon] : sum of water and ice content kg/kg air
        w [lev,lat,lon] : vertical motion m/s
        lh [lev] : approximate level height above ground level (m)
        ff [lat,lon] : fire front field
        three_by_three : True if you want 3x3 subset of default 4x4 plot
    '''
    
    # colour bar stuff
    cmap=plotting._cmaps_['verticalvelocity']
    norm=colors.SymLogNorm(0.25,base=2.)
    cbarform=ticker.ScalarFormatter()
    contours=np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0]))
    
    # get plot extent
    extent = [lon[0],lon[-1],lat[0],lat[-1]]
    plt.close()
    f=plt.figure(figsize=[11,10])
    m_lvls = np.array([10, 20, 32,  40,  48,  56,  64,  72,  80,  88,  96,
       104, 108, 112, 116, 120],int)
    nrows,ncols=4,4
    if three_by_three:
        nrows,ncols=3,3
        m_lvls = m_lvls[np.array([0,3,6,8,10,12,13,14,15])]
    #m_lvls = np.round(np.linspace(0,120,16)).astype(int)
    for i in range(len(m_lvls)):
        plt.subplot(nrows,ncols,i+1)
        
        # top panel is wind speed surface values

        #cs,cb=plotting.map_contourf(extent,w[m_lvls[i]],lat,lon,
        cs,cb=plotting.map_contourf(w[m_lvls[i]],lat,lon,
                                    levels = contours,
                                    cmap=cmap, norm=norm,
                                    cbar=False,)#clabel='m/s', cbarform=cbarform)
        
        plt.title("~ %5.0f m"%lh[m_lvls[i]])
        
        ## Add contour where clouds occur
        #with warnings.catch_warnings():
            # ignore warning when there are no clouds:
            #warnings.simplefilter('ignore')
        plt.contour(lon,lat,qc[m_lvls[i]],np.array([cloud_threshold]), 
                        colors='k')
        
        # add nearby towns
        if extentname is not None:
            plotting.map_add_locations_extent(extentname,hide_text=True)
            
        # Add fire front
        plotting.map_fire(ff,lat,lon,
                          colors='orange',
                          linewidths=3, 
                          alpha=0.5 + 0.5*(i==0),
                          )
    
    # add colour bar
    axes=[0.87, 0.20, 0.04, 0.6] 
    #f.subplots_adjust(top=0.95)
    f.subplots_adjust(wspace=0.001)
    f.subplots_adjust(right=axes[0]-0.03)
    cbar_ax = f.add_axes(axes)
    f.colorbar(cs, cax=cbar_ax, format=cbarform)
    return f


def fig2():
    """
    averaged area time series
    A: plan view showing area averaged
    B: time series of ...   
    C: time series of ...
    """
    weather_series()

def fig3():
    """
    as figure 2, but for day 2
    """
    waroona_day2zoom = [115.65,116.03, -33.1,-32.8]
    single_contour=datetime(2016,1,7,14)
    weather_series(
            day2=True,
            extent=waroona_day2zoom,
            showfirelinehour=single_contour,
            mapname='Waroona_day2.tiff',
    )




def fig4():
    """
    Show isochrones over waroona, look at flux over yarloop and compare to iso
    2 panels: 
        A: isochrones, 
        B: FF contours (hourly coloured by localtime day)
    """
    mr='waroona_run3'
    fig = plt.figure(figsize=[10,10])
    ax1 = plt.subplot(2,1,1)
    ## Top left is the isochrones picture
    iso = image.imread('data/Waroona_Fire_Isochrones.png')
    ax1.imshow(iso)
    plt.xticks([],[])
    plt.yticks([],[])

    # Read fire output
    # area affected by fire
    extentA = [115.65,116.21, -33.05,-32.8]
    
    FFront, SHeat, = fio.read_fire(model_run=mr, dtimes=None, 
                                   extent=extentA, 
                                   firefront=True, 
                                   sensibleheat=True,
                                   day1=True, day2=True
                                   )

    ## PANEL B: firefronts hourly contours
    _,ax2 = fireplan(FFront, 
                     show_cbar=False, 
                     fig=fig,
                     subplot_row_col_n=[2,1,2],
                     color_by_day={'6':'pink','7':'orange','8':None},
                     )
    
    plt.subplots_adjust(left=.05, right=.95, top=.96, bottom=.04,
                        wspace=.03, hspace=.01,
                        )

    
    ## REVIEW: rotating y axis, removing label
    plt.sca(ax2)
    plt.yticks(rotation=77.0)
    plt.xlabel("")
    plt.ylabel("")
    #add_labels(fig)
    fio.save_fig(mr,_sn_,"fig4.png",plt)




def fig5():
    '''
    create vert motion slices
    for 2016 06 14:30 LT (06:30 UTC)
    '''

    dtime=datetime(2016,1,6,6)
    model_run='waroona_run3'
    extentname="waroona"
    HSkip=None
    extent=constants.extents[extentname]
    
    # Read vert motion, clouds
    cubes = fio.read_model_run(model_run, 
                               fdtime=dtime, 
                               extent=extent,
                               HSkip=HSkip)
    w,  = cubes.extract('upward_air_velocity')
    qc, = cubes.extract('qc')
    
    
    ## fire front
    ff_dtimes = utils.dates_from_iris(w)
    ff, = fio.read_fire(model_run=model_run, 
                        dtimes=ff_dtimes, 
                        extent=extent, 
                        firefront=True,
                        HSkip=HSkip)
    
    lh = utils.height_from_iris(w)
    lat = w.coord('latitude').points
    lon = w.coord('longitude').points
    offset=8
    
    #for i in range(len(ff_dtimes)):
    # just want 3rd index (half past hour)
    i=2
        
    utc = ff_dtimes[i]
    ltime = utc+timedelta(hours=offset)
    stitle = ltime.strftime("%Y %b %d %H:%M (LT)")
        
    # 3x3 plot creation    
    fig = vert_motion_slices(qc[i].data, w[i].data,lh,lat,lon,
                           ff=ff[i].data,
                           extentname=extentname,
                           three_by_three=True,
                           )                
    # Save figure
    plt.suptitle(stitle)
    add_labels(fig,skips=[9])
    fio.save_fig(model_run,_sn_,"fig5",plt,)
        


def fig6():
    """
    vertical motion of PCB at 1630 LT
    A) plan view over multiple levels averaged
    B) left to right cross section
    C,D) cross sections diagonally
    """
    model_run='waroona_run3'
    dtime=datetime(2016,1,6,6)
    localtime=True
    HSkip=None
    show_windstream=False
    #print("INFO: starting run_X_transect(",model_run,dtime,")")
    ### First use datetime and extentname to read correct outputs:
    extentname=model_run.split('_')[0]
    extent = constants.extents[extentname]
    
    ## read um output over extent [t, lev, lat, lon]
    cubes = fio.read_model_run(model_run, fdtime=[dtime], extent=extent,
                               add_z=True, add_winds=True, add_topog=True,
                               HSkip=HSkip,)
    
    w, = cubes.extract('upward_air_velocity')
    u, = cubes.extract('u')
    qc, = cubes.extract('qc')
    topog, = cubes.extract('surface_altitude')
    z, = cubes.extract('z_th') 
    
    lat = w.coord('latitude').points
    lon = w.coord('longitude').points

    ## fire front
    ffdtimes = utils.dates_from_iris(w)
    ff,sh = fio.read_fire(model_run=model_run, dtimes=ffdtimes, extent=extent,
                        firefront=True, sensibleheat=True,
                        HSkip=HSkip,)
    
    ## Make transects based on PCB occurrence listed latlons
    # sirivan run transects are twice as wide as waroona (bigger fire)
    lat_pcb, lon_pcb = pcb_occurrences(model_run,ffdtimes)[0]
    X1 = [[lat_pcb,lon_pcb - .14], [lat_pcb, lon_pcb+.14]]
    X2 = [[lat_pcb-0.08, lon_pcb-0.08], [lat_pcb+0.08, lon_pcb+0.08]]
    X3 = [[lat_pcb+0.08, lon_pcb-0.08], [lat_pcb-0.08, lon_pcb+0.08]]
    
    # take mean of vert motion between lvls 25-48 approx 500m - 1500m
    with warnings.catch_warnings():
        # ignore warning from collapsing non-contiguous dimension:
        warnings.simplefilter('ignore')
        wmean = w[:,25:48,:,:].collapsed('model_level_number', iris.analysis.MEAN)

    h0,h1 = utils.height_from_iris(wmean,bounds=True)[0]
    
    # pull data from masked arrays from cubes
    zi = z.data
    topogi = topog.data
    # for each timestep:
    # REVIEW: just look at half past
    for i in [2]:
        ## Pull out data for each time step
        qci = qc[i].data
        wi = w[i].data
        ui = u[i].data
        wmeani = wmean[i].data 
        ffi = None
        shi = None
        if ff is not None:
            ffi = ff[i].data
        if sh is not None:
            shi = sh[i].data
        
        # datetime timestamp for file,title
        labeltime=ffdtimes[i]
        labeltimezone = 'LT' if localtime else 'UTC'
        if localtime:
            offset=fio.run_info[model_run]['UTC_offset']
            labeltime = labeltime+timedelta(hours=offset)
            
        stitle = labeltime.strftime("Vertical motion %%Y %%b %%d %%H:%%M (%s)"%labeltimezone)
        wmeantitle='Mean (%3.0fm - %4.0fm)'%(h0,h1)
        
        plot_X_transect(w=wi, u=ui, qc=qci, z=zi, wmean=wmeani, topog=topogi,
               lat=lat, lon=lon, transect1=X1, transect2=X2, transect3=X3,
               ff=ffi,sh=shi,
               wmeantitle=wmeantitle,
               extentname=extentname,
               show_windstream=show_windstream)

        # Save figure into animation folder with numeric identifier
        plt.suptitle(stitle)
        plt.subplots_adjust(hspace=0.08)
        add_labels(plt.gcf(),skips=[4])
        fio.save_fig(model_run, _sn_, "fig6", plt,
                     dpi=600)


def fig8():
    """
    Jet pull down cross sections, left column is coupled vs right column uncoupled
    A,B) top down with wind streams
    C,D) cross section winds
    E,F) cross section winds
    """
    run1="waroona_run3"
    run2="waroona_run3uc"
    # just looking at 2120 LT, jan 6
    hours=[22,]
    extent=[115.8,116.1,-32.92,-32.82]
    ztop=1000
    columntitles=['coupled','uncoupled']
    subsubdir=None
    
    
    ltoffset=fio.run_info[run1]['UTC_offset']
    if extent is None:
        extent = constants.extents[run1.split('_')[0]]
    nrows=3

    # transects = extent centre +- frac of extent width/height
    y0,x0 = (extent[2]+extent[3])/2.0, (extent[0]+extent[1])/2.0
    # do middle point +- 1/3 of distance to extent edge
    dx = (extent[1]-extent[0]) / 3.0
    dy = (extent[3]-extent[2]) / 3.0
    # transects: [[[lat,lon],[lat,lon]],...]
    transects = [[[y0+dy,x0-dx], [y0+dy, x0+dx]],
                 [[y0,x0-dx], [y0, x0+dx]],
                 ]
    tcolors = ['yellow','blue',]
    

    dtimes=fio.run_info[run1]['filedates'][np.array(hours)]
    
    ## for each hour
    #for dti, dt in enumerate(dtimes):
    # REVIEW: just doing the one time slice
    dti=0
    dt=dtimes[0]

    ## read hour of data for both runs
    cubes1 = fio.read_model_run(run1, fdtime=[dt], extent=extent, 
                                add_topog=True, add_winds=True,
                                add_z=True, add_theta=True)
            
    u1,v1,w1,z1 = cubes1.extract(["u","v","upward_air_velocity","z_th"])
    ctimes=utils.dates_from_iris(w1)
    zd1 = z1.data.data
    theta1, = cubes1.extract("potential_temperature")
    topog1=cubes1.extract("surface_altitude")[0].data
    
    # read fire outputs
    ff1,sh1,u101,v101 = fio.read_fire(model_run=run1,
                                  dtimes=ctimes, 
                                  extent=extent,
                                  sensibleheat=True,
                                  wind=True)
    
    cubes2 = fio.read_model_run(run2, fdtime=[dt], extent=extent, 
                               add_topog=True, add_winds=True,
                               add_z=True, add_theta=True)
    u2,v2,w2,z2 = cubes2.extract(["u","v","upward_air_velocity","z_th"])
    zd2 = z2.data.data
    theta2, = cubes2.extract("potential_temperature")
    topog2=cubes2.extract("surface_altitude")[0].data
    # read fire outputs
    ff2,sh2,u102,v102 = fio.read_fire(model_run=run2,
                                  dtimes=ctimes, 
                                  extent=extent,
                                  sensibleheat=True,
                                  wind=True)
    
    
    
    #for cti,ct in enumerate(ctimes):
    # REVIEW: just doing one time slice
    cti=1
    ct=ctimes[1]
    
    ## FIGURE BEGIN:
    fig = plt.figure(figsize=[9,10])
    ## looking at two separate runs
    for runi, [sh,ff,u10,v10,u,v,w,zd,theta,topog] in enumerate(zip([sh1,sh2],[ff1,ff2],[u101,u102],[v101,v102],[u1,u2],[v1,v2],[w1,w2],[zd1,zd2],[theta1,theta2],[topog1,topog2])):

        # read datacubes
        lats = ff.coord('latitude').points
        lons = ff.coord('longitude').points
        shd = sh[cti].data.data
        LT = ct + timedelta(hours=ltoffset)
        ffd = ff[cti].data.data
        u10d = u10[cti].data.data
        v10d = v10[cti].data.data
        
        #plt.subplot(4,2,1+runi) # top row
        _,ax = topdown_emberstorm(
            fig=fig, subplot_row_col_n=(nrows,2,1+runi),
            extent=extent, lats=lats, lons=lons,
            ff=ffd, sh=shd, u10=u10d, v10=v10d,
            topog=topog,
            annotate=False,
            )
        #ax.set_aspect('equal')
        plt.xticks([115.8,116.1],[115.8,116.1])
        if runi==1:
            plt.yticks([],[])
        else:
            plt.yticks([-32.92,-32.82],[-32.92,-32.82])

        plt.title(columntitles[runi])
        for transect, tcolor in zip(transects,tcolors):
            ## Add dashed line to show where transect will be
            start,end = transect
            # outline the dashed line
            line_effects=[patheffects.Stroke(linewidth=5, foreground='darkgrey'), patheffects.Normal()]
            plt.plot([start[1],end[1]],[start[0],end[0], ], 
                    linestyle='--',
                    color=tcolor,
                    linewidth=3, 
                    path_effects=line_effects,
                    #alpha=0.9,
                    )
        
        ## Plot title
        plt.suptitle(LT.strftime('%b %d, %H:%M (UTC+8)'),fontsize=16)
        
        # Transect plots
        for trani,transect in enumerate(transects):
            ax=plt.subplot(nrows,2,3+runi+2*trani)
            
            trets = transect_emberstorm(
                u[cti].data.data,
                v[cti].data.data,
                w[cti].data.data,
                zd, lats, lons, transect,
                topog=topog,
                sh=shd,
                theta=theta[cti].data.data,
                ztop=ztop,
                theta_contourargs={'levels':[297,300,303,306,309],
                    'linewidths':[3],},
                wind_contourargs={'colorbar':False},
                )
            if runi==1:
                plt.yticks([],[])
            
            # match splines to transect lines
            plotting.set_spine_color(ax,tcolors[trani])
            # finally add desired annotations
            plotting.annotate_max_winds(trets['s'])
            # REVIEW: added distance note
            if trani == 1:
                plt.xlabel("~ 19 km")
            
    # final touches
    plt.tight_layout(rect=[0,0.03,0.92,0.975]) #(left, bottom, right, top)
    ## add wind speed colorbar
    cbar_ax = fig.add_axes([0.91, 0.35, 0.02, 0.3]) # X Y Width Height
    cmap = matplotlib.cm.get_cmap(name=plotting._cmaps_['windspeed'])
    #cmap.set_over('1.0')
    #norm = matplotlib.colors.Normalize(vmin=0, vmax=22.5)
    norm = matplotlib.colors.BoundaryNorm(np.arange(0,22.51,2.5), cmap.N)
    
    cb1 = matplotlib.colorbar.ColorbarBase(cbar_ax, cmap=cmap,
                                           norm=norm,
                                           extend='max',
                                           orientation='vertical')
    cb1.set_label('ms$^{-1}$')
    
        
    ## SAVE FIGURE
    fio.save_fig(run1,_sn_,"fig8",plt=plt,dpi=600)




def fig9(mr="waroona_run3"):
    """
    evening intensity plan views
    A,B) Plan views of intensity
    """
    extent=[115.72,116.03,-33.03,-32.84]
    # just doing topdown emberstorm for 2 time slices:
    t0 = datetime(2016,1,7,11) 
    t0i = 4
    t1 = datetime(2016,1,7,8)
    t1i = 0
    t2 = datetime(2016,1,7,7)
    t2i = 5
    #t3 = datetime(2016,1,7,11)
    #t3i = 4
    wmap_height=300
    for hour,hi in zip([t2,t0,t1],[t2i,t0i,t1i]):
        
        cubes = fio.read_model_run(mr, fdtime=[hour], extent=extent, 
                                   add_topog=True, add_winds=True,
                                   add_z=True, add_theta=True)
        u,v,w,z = cubes.extract(["u","v","upward_air_velocity","z_th"])
        theta, = cubes.extract("potential_temperature")
        topog=cubes.extract("surface_altitude")[0].data
        topogd = topog
        ctimes = utils.dates_from_iris(w)
        
        # extra vert map at ~ 300m altitude
        levh = utils.height_from_iris(w)
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
        # JUST DOING one time slice
        dt=ctimes[hi]
        for dti in [hi]:
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
                
            ## Plot title
            fname=LT.strftime('fig9_%H%M_')+mr
            plt.title(LT.strftime('%b %d, %H:%M (UTC+8)'))
            plt.tight_layout()
            fio.save_fig(mr,_sn_,fname,plt=plt)


if __name__=='__main__':
    
    ### Run the stuff

    ## ~ 35GB ram and 9 min on qsub for fig1 and fig2 together
    ## feedback recieved, some changes to day2 and location font size
    ## just using W and Y for waroona and yarloop!
    #fig2()
    #fig3()

    ## ~ 5GB, 30 seconds on qsub
    ## ready for feedback
    #fig4()

    ## runs without compute node fine < 1 min
    ## ready for feedback
    #fig5()

    ## 
    #fig6()

    fig8()

    #fig9("waroona_run3_day2_early")
    #fig9("waroona_run3")


