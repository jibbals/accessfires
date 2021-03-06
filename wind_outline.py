#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 20:45:05 2019
    Winds outline script
        Transect views of horizontal and vertical winds
@author: jesse
"""

import matplotlib
# don't plot on screen, send straight to file
# this is for lack of NCI display
matplotlib.use('Agg',warn=False)

# plotting stuff
from matplotlib import colors, ticker
import matplotlib.pyplot as plt
#import matplotlib.ticker as tick
#import matplotlib.patches as mpatches
import numpy as np
import cartopy.crs as ccrs # projection stuff
from datetime import timedelta
import warnings

# local modules
from utilities import plotting, utils, fio, constants

###
## GLOBALS
###
_sn_ = 'wind_outline'

_transects_             = {} 
__x0__,__x1__ = 115.8, 116.19

_transects_['waroona1'] = [-32.79   , __x0__], [-32.92   , __x1__]
_transects_['waroona2'] = [-32.82   , __x0__], [-32.93   , __x1__]
_transects_['waroona3'] = [-32.86   , __x0__], [-32.88   , __x1__]
# again but start more southerly and end more northerly
_transects_['waroona4'] = [-32.92   , __x0__], [-32.82   , __x1__]
_transects_['waroona5'] = [-32.96   , __x0__], [-32.85   , __x1__] 
_transects_['waroona6'] = [-32.87   , __x0__], [-32.89   , __x1__]

# looking at sir ivan
#_extents_['sirivan']    = [149.2, 150.4, -32.4, -31.6]
_transects_['sirivan1'] = [-32.05, 149.4  ], [-32.0  , 150.3 ]
_transects_['sirivan2'] = [-32.0 , 149.4  ], [-31.95 , 150.3 ]
_transects_['sirivan3'] = [-32.1 , 149.4  ], [-32.15 , 150.3 ]
_transects_['sirivan4'] = [-31.8 , 149.45 ], [-32.15 , 150.2 ]
_transects_['sirivan5'] = [-31.95, 149.4  ], [-31.80 , 150.2 ]
_transects_['sirivan6'] = [-31.7 , 149.5  ], [-32.1  , 150.3 ]


def show_transects():
    """
    For sirivan and waroona show the transects on a contour map with final fire outline as well
    """
    
    # Show transects
    waroona_xticks = np.arange(115.8, 116.21, .1)
    waroona_yticks = np.arange(-33, -32.701, .05)
    sirivan_xticks = np.arange(149.2, 150.4, .2) #149.2, 150.4
    sirivan_yticks = np.arange(-32.4,-31.6,.1) #-32.4, -31.6
    for extentname, xticks, yticks in zip(['waroona','sirivan'],
                                          [waroona_xticks,sirivan_xticks],
                                          [waroona_yticks,sirivan_yticks]):
        extent = constants.extents[extentname]
        mr = '%s_run1'%extentname
        
        # read topog, fire
        topog = fio.read_topog(mr,extent=extent)
        lat,lon = topog.coord('latitude').points, topog.coord('longitude').points
        ff, = fio.read_fire(mr,extent=extent,
                            dtimes=[fio.run_info[mr]['filedates'][-1]])

        plotting.map_topography(extent,topog.data,lat,lon)
        plt.title("Transects")
        ## show fire contour if model run has fire
        if ff is not None:
            plt.contour(lon,lat, np.transpose(ff[0].data),np.array([0]),colors='red')
        
        ## plot each transect 
        for transect in range(6):
            start,end = _transects_["%s%d"%(extentname,transect+1)]
            plt.plot([start[1],end[1]],[start[0],end[0], ], '--',
                     linewidth=2, label='X%d'%(transect+1),
                     marker='X', markersize=7,markerfacecolor='white')
        
        # add nearby towns
        plotting.map_add_locations_extent(extentname)
        
        plt.legend()
        plt.yticks(yticks)
        plt.xticks(xticks)
        pname="figures/transects_%s.png"%extentname
        fio.save_fig_to_path(pname,dpi=600,plt=plt)
    
    
def wind_outline(s,u,v,w,
                 qc, topog,
                 z,lat,lon,
                 transect,
                 ff=None,
                 extentname='waroona',
                 ztop = 5000,
                 cloud_threshold=constants.cloud_threshold,
                ):
    '''
    311 Plot showing contour map and wind speed, along with near sites and transect
    312 plot showing vertical motion along transect
    313 plot showing wind speed along transect
    INPUTS:
        wind speed, x-wind, y-wind, vert-wind, cloud content, topography, model level altitude, lat,lon, datetime (for titles)
        ff: fire front array, negative where fire has burned, optional
        extentname = { 'waroona' | 'sirivan' } choice of two fire locations
        transect = [[lat0,lon0],[lat1,lon1]]
        nquivers: int, roughly how many quivers are desired
        quiverscale: int, changes how long the arrows are (also may need to fiddle)
        cloud_threshold: float, g/kg where cloud outline is drawn
    '''
    # set font sizes
    plotting.init_plots()
    
    # get plot extent, and transect
    extent = [lon[0],lon[-1],lat[0],lat[-1]]
    start,end = transect
    
    plt.figure(figsize=[7,10])
    plt.subplot(3,1,1)
    
    # top panel is topography
    plotting.map_topography(extent,topog,lat,lon)
    plt.title('Topography, winds')
    
    # Add fire front contour
    if ff is not None:
        plotting.map_fire(ff,lat,lon, transform=False)
    
    # start to end x=[lon0,lon1], y=[lat0, lat1]
    plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
             linewidth=2, 
             marker='>', markersize=7,markerfacecolor='white')
    
    # add nearby towns
    if extentname is not None:
        plotting.map_add_locations_extent(extentname)
    
    ## add streamplot
    plt.get_cmap(plotting._cmaps_['windspeed'])
    plt.streamplot(lon,lat,u[0],v[0],color='darkslategrey',)
    
    ## Second row is transect plots
    plt.subplot(3,1,2)
    nx=utils.number_of_interp_points(lat,lon,start,end)
    wslice, slicex, slicez  = plotting.transect_w(w,z, lat, lon,start,end,
                                                  topog=topog,
                                                  ff=ff,
                                                  npoints=nx,
                                                  ztop=ztop,
                                                  lines=None)
    plt.ylabel('height (m)')
    ## Add contour where clouds occur
    qcslice, _,_ = utils.transect(qc,lat,lon,start,end,nx=nx)
    with warnings.catch_warnings():
        # ignore warning when there are no clouds:
        warnings.simplefilter('ignore')
        plt.contour(slicex,slicez,qcslice,np.array([cloud_threshold]),colors='teal')
    
    ax3 = plt.subplot(3,1,3)
    trs, trx, trz = plotting.transect_s(s,z,lat,lon,start,end,
                                        topog=topog, ztop=ztop,
                                        ff=ff, npoints=nx)
    xticks,xlabels = plotting.transect_ticks_labels(start,end)
    plt.xticks(xticks[0::2],xlabels[0::2]) # show transect start and end
    
    # Annotate max wind speed
    # only care about lower troposphere
    # 60 levels is about 2800m, 80 levels is about 4700m, 70 levels : 3700m
    upto = 70 
    mloc = np.unravel_index(np.argmax(trs[:upto,:],axis=None),trs[:upto,:].shape)
    note="max = %5.1f m/s\n  (at %4.0f m) "%(trs[:upto,:][mloc], trz[:upto,:][mloc])
    trans = ax3.get_xaxis_transform() # x in data untis, y in axes fraction
    ax3.annotate(note, fontsize=15,
                 xy=(0.33, -0.2 ), xycoords=trans)
    
    return plt

def outline_model_winds(model_run='sirivan_run1', hours=None, dpi=200):
    """
    Run wind_outline for one output file from a chosen 'model_run'
    """
    
    extentname=model_run.split('_')[0]
    extent=constants.extents[extentname]
    
    ## Read the cubes
    cubes = fio.read_model_run(model_run, fdtime=hours, extent=extent,
                               add_topog=True,
                               add_z=True, 
                               add_winds=True)
    
    zth, topog = cubes.extract(['z_th','surface_altitude'])
    u,v,s, qc, w = cubes.extract(['u','v','s','qc','upward_air_velocity'])
    cubetimes = utils.dates_from_iris(u)
    
    ## fire front
    ff, = fio.read_fire(model_run=model_run, dtimes=cubetimes, extent=extent, 
                        firefront=True)
    lat, lon = w.coord('latitude').points, w.coord('longitude').points
    
    ## loop over available timedate steps
    for i, dt in enumerate(cubetimes):
        stitle = dt.strftime("%Y %b %d %H:%M (UTC)")
        si, ui, vi, wi, qci, topogi, zthi, ffi = [s[i].data.data, u[i].data.data, 
                                                 v[i].data.data, w[i].data.data, 
                                                 qc[i].data.data, topog.data.data, 
                                                 zth.data.data, None]
        if ff is not None:
            ffi=ff[i].data.data

        ## loop over different transects
        for ii in range(6):
            transect = _transects_["%s%d"%(extentname,ii+1)]
            plt = wind_outline(si, ui, vi, wi, qci, topogi, zthi, lat, lon,
                               transect=transect,
                               ff = ffi,
                               extentname=extentname)
            # Save figure into subfolder with transect identifier
            plt.suptitle(stitle)
            fio.save_fig(model_run, _sn_, cubetimes[i], plt,
                         subdir='transect_%d'%(ii+1),
                         dpi=dpi)

def vorticity_layers(model_run="waroona_run2", hour=16, levels=[3,5,10,20,30,40],
                     vorticity_max=0.01,
                     extent=None, HSkip=None):
    '''
    Show horizontal slices of vorticity and horizontal winds on several vertical levels. 
    '''
    extentname=None
    if extent is None:
        extentname=model_run.split('_')[0]
        extent = constants.extents[extentname]
        
    fdtime = fio.run_info[model_run]['filedates'][hour]
    
    # read cubes
    cubes = fio.read_model_run(model_run, fdtime=[fdtime], extent=extent, 
                               add_winds=True,
                               HSkip=HSkip)
    u,v = cubes.extract(['u','v'])
    w, = cubes.extract('upward_air_velocity')
    lat = w.coord('latitude').points
    lon = w.coord('longitude').points
    height = utils.level_from_height(w)
        
    dtimes = utils.dates_from_iris(u)
        
    ff, = fio.read_fire(model_run, dtimes, extent=extent, HSkip=HSkip)
    
    # constant colorbar for vorticity
    cmap1="PuOr"
    cmap2="viridis"
    cmap3="viridis"
    contours1=np.linspace(-1*vorticity_max,vorticity_max,50)
    contours2=np.linspace(0,1,50)
    contours3=np.linspace(0,10,50)
    
    # streamplot density of lines
    density_x,density_y = .5,.3
    
    for di, dtime in enumerate(dtimes):
        fig = plt.figure(figsize=[10,10])
        
        for ii,lev in enumerate(levels):
            
        
            Ui = u[di,lev].data.data
            Vi = v[di,lev].data.data
            zi,OWi,OWZi = utils.vorticity(Ui,Vi,lat,lon)
            
            # plot vorticity (units will be wrong: dims are degrees not metres)
            for jj, metric in enumerate([zi,OWi,OWZi]):
                plt.subplot(len(levels),3,ii*3+jj+1)
                
                #metric[~np.isfinite(metric)] = np.NaN
                
                if jj == 0:
                    cs1 = plt.contourf(lon,lat, metric, contours1, cmap=cmap1, extend='both')
                elif jj == 1:
                    print("DEBUG1:",metric.shape,np.nanmin(metric),np.nanmax(metric))
                    metric[metric<0]=np.NaN
                    cs2 = plt.contourf(lon,lat, metric, contours2, cmap=cmap2, extend='max')
                elif jj == 2:
                    print("DEBUG2:",metric.shape,np.nanmin(metric),np.nanmax(metric))
                    metric[metric<0]=np.NaN
                    cs3 = plt.contourf(lon,lat, metric, contours3, cmap=cmap3, extend='max')
            
                # show winds
                plt.streamplot(lon,lat,Ui,Vi, color='k', 
                               density=(density_x, density_y))
            
                if extentname is not None:
                    plotting.map_add_locations_extent(extentname, hide_text=True)
        
                plt.xticks([],[])
                plt.yticks([],[])
                if jj == 0:
                    plt.ylabel("%.0f m"%(height[lev]),fontsize=13)
            
            
                # add fire contour
                if ff is not None:
                    plotting.map_fire(ff[di].data.data,lat,lon)
        
                if ii==0:
                    plt.title(['vorticity (1/s)','OW (1)','OWZ (s)'][jj])
        
        # reduce vert gap between subplots
        fig.subplots_adjust(hspace=0.1)
        
        # add colourbar for vorticity:
        cbar_ax = fig.add_axes([0.1, 0.05, 0.2, 0.05])# X Y Width Height
        cb1 = fig.colorbar(cs1, cax=cbar_ax, orientation='horizontal',
                           format=matplotlib.ticker.ScalarFormatter(), 
                           pad=0)
        
        # add colourbar for normalised metrics:
        cbar_ax2 = fig.add_axes([0.4, 0.05, 0.2, 0.05])# X Y Width Height
        cb2 = fig.colorbar(cs2, cax=cbar_ax2, orientation='horizontal',
                           format=matplotlib.ticker.ScalarFormatter(), 
                           pad=0)
        
        # add colourbar for vorticity:
        cbar_ax3 = fig.add_axes([0.7, 0.05, 0.2, 0.05])# X Y Width Height
        cb3 = fig.colorbar(cs3, cax=cbar_ax3, orientation='horizontal',
                           format=matplotlib.ticker.ScalarFormatter(), 
                           pad=0)
        
        tv1 = [-vorticity_max, 0, vorticity_max]
        tv2 = [0,0.5,1]
        tv3 = [0,5,10]
        for cb,tv in zip([cb1,cb2,cb3],[tv1,tv2,tv3]):
            cb.set_ticks(list(tv))
            cb.ax.set_xticklabels(list(tv))
        
        # save figure
        fio.save_fig(model_run=model_run, script_name=_sn_,pname=dtime,
                     plt=plt,subdir="vorticity")

def transects_hwinds(model_run, hour=18, transects=None, extent=None, ztop=4000,
                     subdir="hwinds", HSkip=None):
    """
    Arguments:
        hour: int from 0 to 23 (or 47 if second day is available)
        transects: [[lat,lon0,lon1], ...]
        extent: can change mapping extent if desired
        ztop: height of top of transect
        subdir: subfolder for plots
        HSkip: reduce horizontal resolution for RAM saving
    Figure:
        Region, with transect lines, firefront, and 10m winds
        row for each transect showing horizontal winds 
    """
    extentname=model_run.split('_')[0]
    if extent is None:
        extent=constants.extents[extentname]
    nrows=len(transects)+1
    model_hours = fio.run_info[model_run]['filedates']
    dt = model_hours[hour]
    # First read data
    cubes = fio.read_model_run(model_run, fdtime=dt,
                               HSkip=HSkip, 
                               add_winds=True, add_z=True, add_topog=True)
    u,v,s,z,topog,w = cubes.extract(['u','v','s','z_th','surface_altitude','upward_air_velocity'])
    ctimes = utils.dates_from_iris(u)
    clats, clons = u.coord('latitude').points, u.coord('longitude').points
    
    ff,u10,v10 = fio.read_fire(model_run, dtimes=ctimes, extent=extent, HSkip=HSkip, 
                               wind=True)
    flats,flons = ff.coord('latitude').points, ff.coord('longitude').points
    
    # loop over time slices
    for di, dtime in enumerate(ctimes):
        # figure:
        fig = plt.figure(figsize=[12,12])
        _, ax1 = plotting.map_tiff_qgis(fname=extentname+'.tiff', 
                                         extent=extent, 
                                         fig=fig, 
                                         subplot_row_col_n=[nrows,1,1],)
        
        # plot firefront
        plotting.map_fire(ff[di].data,flats,flons)
        # add winds streamplot
        plt.streamplot(flons,flats,u10[di].data,v10[di].data, color='white',
                       density=(0.8, 0.5),
                       )
        plt.title('10m horizontal winds')
        # loop over transects
        for ti, [lat,lon0,lon1] in enumerate(transects):
            # add transect to map
            start,end = [lat,lon0],[lat,lon1]
            plt.sca(ax1)
            plt.plot([lon0,lon1],[lat,lat], '--', color='k', 
                     linewidth=2, 
                     ) 
            
            # plot transect of horizontal wind speed on new row
            plt.subplot(nrows,1,ti+2)
            npoints=utils.number_of_interp_points(clats,clons,start,end)
            # contourf of horizontal wind speeds
            _, slicex, slicez = plotting.transect_s(s[di].data, z.data, 
                                                    clats,clons, 
                                                    start, end,
                                                    lines=None, npoints=npoints,
                                                    topog=topog.data)
            
            # east-west and vertical winds on transect
            uslice,_,_ = utils.transect(u[di].data,clats,clons,start,end,nx=npoints)
            wslice,_,_ = utils.transect(w[di].data,clats,clons,start,end,nx=npoints)
            
            # Streamplot
            plotting.streamplot_regridded(slicex,slicez,uslice,wslice,
                                          density=(1.5,1.5), 
                                          color='darkslategrey',
                                          zorder=1)
            
            plt.title("lat: %.3f"%lat)
        plt.suptitle(dtime.strftime("Horizontal winds %Y%m%d %H%M (UTC)"))
        fio.save_fig(model_run, script_name=_sn_, pname=dtime, plt=plt, 
                     subdir=subdir)
 

def vertical_vortex(mr='waroona_run3',
                    zlevels=[5,10,15,25,35,45],
                    hours=[12,13,14],
                    extent=None,
                    HSkip=None,
                    minimap=True,
                    **map_tiff_args
                    ):
    """
    show horizontal slices of horizontal wind streamplots 
    with vertical motion contourf and optional minimap
    layout:
        [ 3 ] [ 6 ]
        [ 2 ] [ 5 ]
        [ 1 ] [ 4 ]
        [ minimap ]
    """
    ##### set up extent, colormap, streamplot linewidths
    extentname=mr.split('_')[0]
    if extent is None:
        extent = constants.extents[extentname]
    dtimes = fio.run_info[mr]['filedates'][hours]
    # vertical wind colourbar is constant
    wcmap=plotting._cmaps_['verticalvelocity']
    wnorm=colors.SymLogNorm(0.25,base=2.) # linear to +- 0.25, then log scale
    wcontours=np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0]))
    speedmax=20
    
    ##### loop over desired hours #####
    for dtime in dtimes:
        
        ##### Grab hour of model run output #####
        cubes = fio.read_model_run(mr, fdtime=dtime,
                                   extent=extent,
                                   HSkip=HSkip, 
                                   add_winds=True, 
                                   add_z=True, 
                                   add_topog=True)
        u,v,s,z,topog,w = cubes.extract(['u','v','s','z_th','surface_altitude','upward_air_velocity'])
        ctimes = utils.dates_from_iris(u)
        lats, lons = u.coord('latitude').points, u.coord('longitude').points
        height = utils.level_from_height(w)
        
        # maybe want minimap with fire etc.
        if minimap:
            # extent plus a little for clarity
            dx=extent[1]-extent[0]
            dy=extent[3]-extent[2]
            full_extent= extent[0]-dx/4.0, extent[1]+dx/4, extent[2]-dy/5.0, extent[3]+dy/5.0
            FF, SH, u10,v10 = fio.read_fire(model_run=mr, 
                                            dtimes=ctimes, 
                                            extent=full_extent, 
                                            firefront=True, 
                                            sensibleheat=True,
                                            wind=True,
                                            )
            flons,flats = FF.coord('longitude').points, FF.coord('latitude').points
        
        ##### Loop over time slices #####
        for cti, ctime in enumerate(ctimes):
            
            fig=plt.figure(figsize=(12,16))
            levelmap = {0:5,1:3,2:1,3:6,4:4,5:2} # map level index to subplot index
            for i in range(6):
                ax=plt.subplot(3+minimap,2,levelmap[i])
                zi=zlevels[i]
                Ui = u[cti,zi].data.data
                Vi = v[cti,zi].data.data
                Si = np.hypot(Ui,Vi)
                Wi = w[cti,zi].data.data
                # Plot vertical motion colours
                cs=plt.contourf(lons, lats, Wi, 
                                levels=wcontours,
                                cmap=wcmap, norm=wnorm)
                # streamplot of horizontal winds
                streamLW=utils.wind_speed_to_linewidth(Si,speedmax=speedmax,lwmax=6)
                plt.streamplot(lons,lats,Ui,Vi,
                               linewidth=streamLW,
                               color='k',
                               density=[.5,.5],
                               )
                # reset axes
                plt.xlim(lons[0],lons[-1])
                plt.ylim(lats[0],lats[-1])
                plt.title("~%.2fkm"%(height[zi]/1000.))
                
                # fix ticks
                if (levelmap[i]%2) == 1:
                    ax.yaxis.set_major_locator(ticker.LinearLocator(numticks=4))
                    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
                else:
                    plt.yticks([],[])
                if levelmap[i] > 4:
                    ax.xaxis.set_major_locator(ticker.LinearLocator(numticks=4))
                    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
                else:
                    plt.xticks([],[])
            
            if minimap:
                
                # defaults:
                if 'fname' not in map_tiff_args:
                    map_tiff_args['fname']=extentname+'.tiff'
                if 'extent' not in map_tiff_args:
                    map_tiff_args['extent']=full_extent
                map_tiff_args['subplot_row_col_n']=[4,1,4]
                map_tiff_args['fig'] = fig
                # grab tiff and plot onto row 1 of 4
                fig,ax = plotting.map_tiff_qgis(**map_tiff_args)
                
                # add heat flux and fire and winds
                plotting.map_sensibleheat(SH[cti].data,flats,flons,alpha=0.6,
                                          colorbar=False, cbar_kwargs={'label':"Wm$^{-2}$"})
                plotting.map_fire(FF[cti].data,flats,flons)
                S10 = np.hypot(u10[cti].data,v10[cti].data)
                streamLW=utils.wind_speed_to_linewidth(S10,speedmax=speedmax,lwmax=5)
                
                plt.streamplot(flons,flats,u10[cti].data,v10[cti].data,
                               linewidth=streamLW,
                               density=[0.4,0.4])
                # reset axes
                plt.xlim(flons[0],flons[-1])
                plt.ylim(flats[0],flats[-1])
                # add rectangle
                plotting.map_add_rectangle(extent,edgecolor='darkgrey')
                plotting.map_add_locations_extent(extentname,hide_text=False)
                ax.yaxis.set_major_locator(ticker.LinearLocator(numticks=5))
                ax.xaxis.set_major_locator(ticker.LinearLocator(numticks=5))
                ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
                ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
            
            offsethours=8 if np.min(lons)<120 else 10
            ltime=ctime+timedelta(hours=offsethours)
            plt.suptitle(ltime.strftime("wind slices at %Y %b %d %H:%M (LT)"))
            
            # add vert wind colourbar
            cbar_ax = fig.add_axes([0.49, 0.4+0.1*minimap, 0.01, 0.2]) # X Y Width Height
            fig.colorbar(cs, cax=cbar_ax, format=ticker.ScalarFormatter(), pad=0)
            plt.tight_layout(rect=[0, 0.03, 1, 0.95]) #LEFT,BOT,Right,Top
            fio.save_fig(mr,_sn_,ctime,plt=plt,subdir="vortex")
        
    
#########################################################################
#########################################################################
#########################################################################


if __name__ == '__main__':
    
    if False:
        #         [lat, lon0, lon1],
        transects = [ [-32.84, 115.82, 116.05], 
                     [-32.9, 115.82, 116.05], 
                     [-32.96, 115.82, 116.05], ]
        extent=[115.6,116.21, -33.05,-32.75]
        for hour in range(14,15):
            transects_hwinds(model_run='waroona_run1', extent=extent,
                             hour=hour,
                             transects=transects)
    
    if False:
        extent=None
        SI_PCB = [149.5,150.2,-32.15,-31.95] # lon0,lon1,lat0,lat1
        Waroona = constants.extents['waroona']
        for hour in np.arange(16,18):
            vorticity_layers("waroona_run1",
                             hour=hour,
                             extent=Waroona,)
    
    if False:
        allmr = fio.run_info.keys()
        allmr = ['waroona_run1']
        for mr in allmr:
            hours = fio.run_info[mr]['filedates'][14:16]
            for hour in hours:
                print("info: wind_outline", mr, hour)
                outline_model_winds(mr, hours=[hour])
    
    if True:
        vertical_vortex(
            mr='sirivan_run1',
            extent=constants.extents['sirivanz'],
            hours=range(5,24))
