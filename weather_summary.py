#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 08:45:14 2019
    Weather summary looking at winds and clouds
@author: jesse
"""
import matplotlib
matplotlib.use("Agg",warn=False)

from matplotlib import colors, ticker, patches
import matplotlib.dates as mdates

import numpy as np
import matplotlib.pyplot as plt
import iris # minimal fio done here 
import cartopy.crs as ccrs
import warnings
from datetime import datetime, timedelta

from utilities import utils, plotting, fio, constants

## GLOBAL
#Script name
_sn_ = 'weather_summary'
__cloud_thresh__ = constants.cloud_threshold

def plot_weather_summary(U,V,W, height, lat, lon, extentname, 
                         Q=None, FF=None, Streamplot=True):
    '''
    Show horizontal slices of horizontal and vertical winds averaged between 
    several vertical levels. Also shows clouds (Q) and fires (FF) with contour outline.
    
    INPUTS: 
        U,V,W: wind speed in lon, lat, vert dimension m/s [z,lat,lon]
        height: altitude ASL or AGL m [z,lat,lon] or [z]
        height: level heights for titles (array [ levs ])
    '''
    
    row1 = (100<=height) * (height<500)
    row2 = (500<=height) * (height<1500)
    row3 = (1500<=height) * (height<3000)
    row4 = (3000<=height) * (height<5000)
    # todo row 5
    row5 = (5000<height) * (height<10000)
    
    # vertical wind colourbar is constant
    wcmap=plotting._cmaps_['verticalvelocity']
    wnorm=colors.SymLogNorm(0.25) # linear to +- 0.25, then log scale
    wcontours=np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0]))
    
    fig = plt.figure(figsize=[10,10])
    for ii,row in enumerate([row1, row2, row3, row4, row5]):
        plt.subplot(5,2,ii*2+1)
        
        Ui = U[row]
        Vi = V[row]
        Si = np.hypot(Ui,Vi)
        
        # mean of levels # .data is masked array, just use base array
        Ur = np.mean(Ui, axis=0)
        Vr = np.mean(Vi, axis=0)
        Sr = np.mean(Si, axis=0)
        
        # wind speed contourf
        plt.contourf(lon, lat, Sr, 10)
        plt.colorbar(ticklocation=ticker.MaxNLocator(5),pad=0)
        
        if extentname is not None:
            plotting.map_add_locations_extent(extentname, hide_text=True)
        
        #Streamplot the horizontal winds
        # TODO: convert lats,lons to metres?
        plt.streamplot(lon,lat,Ur,Vr, 
                       color='k', 
                       density=(0.6, 0.5))
        
        # remove x and y ticks
        plt.xticks([],[])
        plt.yticks([],[])
        # add ylabel
        plt.ylabel("%.0f - %.0f m"%(height[row][0], height[row][-1]),fontsize=13)
        if ii==0:
            plt.title('horizontal winds (m/s)')
        
        # Now plot vertical motion
        plt.subplot(5,2,ii*2+2)
        Wi = W[row]
        Wr = np.mean(Wi,axis=0)
        
        cs = plt.contourf(lon, lat, Wr, wcontours, cmap=wcmap, norm=wnorm)
        if extentname is not None:
            plotting.map_add_locations_extent(extentname, hide_text=ii>0)
        
        # add cloud hatching
        if Q is not None:
            Qmax = np.max(Q[row], axis=0)
            # add hatches over where Qmax is greater than cloud threshhold
            plt.contourf(lon,lat,Qmax, [0,__cloud_thresh__,100], 
                         colors=['None']*3, hatches=[None,'/','//'],
                         extend='both',)
        
        if FF is not None:
            plotting.map_fire(FF, lat, lon)
        
        plt.xticks([],[])
        plt.yticks([],[])
        if ii==0:
            plt.title('vertical motion (m/s)')
        
    # reduce vert gap between subplots
    fig.subplots_adjust(hspace=0.1)
    # add vert wind colourbar
    cbar_ax = fig.add_axes([0.905, 0.4, 0.01, 0.2])# X Y Width Height
    fig.colorbar(cs, cax=cbar_ax, format=ticker.ScalarFormatter(), pad=0)
        

def weather_summary_model(model_version='waroona_run1',
                          fdtimes=None,
                          zoom_in=None,
                          HSkip=None):
    '''
    Read model run output hour by hour, running plot_weather_summary on each
    time slice. Can subset to just show fdtimes, and can also zoom to specific 
    lat lon box, which saves into the 'zoomed' subfolder
    '''
    # font sizes etc
    plotting.init_plots()
    
    extentname = model_version.split('_')[0]
    extent = plotting._extents_[extentname]
    if zoom_in is not None:
        extentname=None
        extent = zoom_in
    if fdtimes is None:
        fdtimes = fio.model_outputs[model_version]['filedates']
    FF = None
    
    # read one hour at a time, plot each available time slice
    for fdtime in fdtimes:
        
        # read cubes
        cubes = fio.read_model_run(model_version, fdtime=[fdtime], extent=extent, 
                                   add_winds=True,
                                   HSkip=HSkip)
        u,v = cubes.extract(['u','v'])
        w, = cubes.extract('upward_air_velocity')
        clouds = cubes.extract('qc')
        qc = None
        if len(clouds) == 1:
            qc = clouds[0]
        lat = w.coord('latitude').points
        lon = w.coord('longitude').points
        height = w.coord('level_height').points
        dtimes = utils.dates_from_iris(u)
        # read fire front
        ff, = fio.read_fire(model_version, dtimes, extent=extent, HSkip=HSkip)
        
        # for each time slice create a weather summary plot
        for i,dtime in enumerate(dtimes):
            
            ui, vi = u[i].data.data, v[i].data.data
            wi = w[i].data.data
            qci = None
            if qc is not None:
                qci = qc[i].data.data
            if ff is not None:
                FF = ff[i].data.data
            
            plot_weather_summary(ui, vi, wi, height, lat, lon, 
                                 extentname=extentname,
                                 Q = qci, FF=FF)
            
            plt.suptitle("%s weather "%model_version + dtime.strftime("%Y %b %d %H:%M (UTC)"))
            subdir=None
            if zoom_in is not None: 
                subdir = 'zoomed'
            fio.save_fig(model_version,_sn_, dtime, plt, subdir=subdir)

def weather_series(model_run='waroona_run1', 
                   extent=None,
                   HSkip=None,):
    """
        time series of surface weather
    """
    
    
    extentname=model_run.split('_')[0]
    # zoomed extent for analysis
    extentnamez = extentname + 'z'
    if extent is None:
        extent=plotting._extents_[extentnamez]
    localtime = False
    proj_latlon = ccrs.PlateCarree() # lat,lon projection
    qc_thresh = constants.cloud_threshold
    
    ## read all the fire data, and lats/lons
    ff,sh, = fio.read_fire(model_run=model_run, dtimes=None,
                           extent=extent, HSkip=HSkip,
                           firefront=True, sensibleheat=True)
    ff=ff[-1] # just final firefront is needed
    
    
    lons,lats = sh.coord('longitude').points, sh.coord('latitude').points
    assert sh is not None, "Missing sensible heat file"
    # datetimes list, convert to local time if desired
    ftimes = utils.dates_from_iris(sh)
    if localtime:
        offset = timedelta(hours=8)
        if extent[0] > 130:
            offset = timedelta(hours=10)
        ftimes = np.array([ft + offset for ft in ftimes ])
    
    # spatial sum over time of sensible heat gives firepower
    # gives GWatts
    firepower = np.sum(utils.firepower_from_cube(sh), axis=(1,2))
    
    # remove zeros:
    prefire = np.isclose(np.cumsum(firepower), 0)
    firepower[prefire] = np.NaN
    
            
    ## Read PFT
    # calculated using kevin's code, stored as GWatts
    lat,lon = plotting._latlons_["fire_%s_upwind"%extentname]
    pft, ptimes, _, _ = fio.read_pft(model_run,lats=lat,lons=lon)
    #pft = np.nanmean(pft,axis=(1,2)) # average spatially
    if localtime:
        ptimes = np.array([pt + offset for pt in ptimes ])
    
    ## Read Temperature:
    # just want surface for most of these, lets just take bottom two rows:
    # surf_constraint = iris.Constraint(model_level_number = lambda cell: cell<=2)
    #cubes=fio.read_model_run(model_run, extent=extent, constraints=surf_constraint, 
    #                         add_RH=True)
    
    # READ EVERYTHING, SUBSET, CALC DESIRED METRICS, PLOT METRICS
    cubes = fio.read_model_run(model_run, extent=extent,
                               HSkip=HSkip)
    ctimes = utils.dates_from_iris(cubes[0])
    
    ## get temperature, RH, cloud
    q,T,qc = cubes.extract(['specific_humidity','air_temperature','qc'])
    
    index_500m = np.argmin(np.abs(q.coord('level_height').points - 500))
    assert index_500m > 0, "Index500 didn't work = %d"%index_500m
    
    clats, clons = q.coord('latitude').points, q.coord('longitude').points
    # compute RH from specific and T in kelvin
    T.convert_units('K')
    # just want surface for Temp and RH
    q = q[:,0,:,:]
    T = T[:,0,:,:].data.data
    qc_sum = qc.collapsed('model_level_number', iris.analysis.SUM)
    qc_frac = np.sum(qc_sum.data > qc_thresh, axis=(1,2))/(len(clats)*len(clons))
    qc_weight = qc_sum.collapsed(['longitude','latitude'], iris.analysis.SUM).data
    #qc_weight[qc_weight <= 0.0001] = np.NaN # take out zeros
    qc_q3 = np.nanpercentile(qc_weight, 75)
    qc_q2 = np.nanpercentile(qc_weight, 50)
    qc_heavy = qc_weight > qc_q3
    qc_mid = (qc_weight > qc_q2) * (qc_weight < qc_q3)
    #qc_weight[np.isnan(qc_weight)] = 0.0 # return the zeros
    RH = utils.relative_humidity_from_specific(q.data.data, T)
    RH = np.mean(RH, axis=(1,2))
    T = np.mean(T, axis=(1,2))
    
    ## get wind speed/ wind dir
    u, v = cubes.extract(['x_wind','y_wind'])
    ## 10m u and v winds (already destaggered) are in fire model output
    u10, v10 = fio.read_fire(model_run=model_run, extent=extent, 
                             dtimes=ctimes,
                             firefront=False, wind=True, 
                             HSkip=HSkip)
    ## wind speed at 10m is output by the fire model
    u10=np.swapaxes(u10.data, 1,2) # time, lon, lat -> time, lat, lon
    v10=np.swapaxes(v10.data, 1,2) # also for v10
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
    ws500_all = utils.wind_speed_from_uv_cubes(u500,v500).data.data
    ws500_q1 = np.quantile(ws500_all, 0.25, axis=(1,2))
    ws500_q3 = np.quantile(ws500_all, 0.75, axis=(1,2))
    #ws500_max = np.max(ws500_all, axis=(1,2))
    ws500 = np.mean(ws500_all, axis=(1,2))
    
    # mean wind direction based on mean u,v
    u10_mean = np.mean(u10,axis=(1,2))
    v10_mean = np.mean(v10,axis=(1,2))
    u500_mean = np.mean(u500.data.data,axis=(1,2)) 
    v500_mean = np.mean(v500.data.data,axis=(1,2))
    #wd0 = utils.wind_dir_from_uv(u10_mean,v10_mean)
    
    #######################
    #### PLOTTING PART ####
    #######################
    
    # figure and multiple y axes for time series
    fig = plt.figure(figsize=[12,12])
    dx,dy = extent[1]-extent[0], extent[3]-extent[2]
    extentplus = np.array(extent)+np.array([-0.3*dx,0.3*dx,-0.2*dy,0.2*dy])
    
    # map with extent shown
    fig, ax_map, proj = plotting.map_tiff_qgis(fname="%s.tiff"%extentname, 
                                               fig=fig,
                                               extent=list(extentplus), 
                                               subplot_row_col_n = [3,1,1],
                                               locnames=[extentname,]
                                               )
    # Add rectangle
    botleft = extent[0],extent[2]
    ax_map.add_patch(patches.Rectangle(xy=botleft,
                                       width=dx,
                                       height=dy,
                                       #facecolor=None,
                                       fill=False,
                                       edgecolor='blue',
                                       linewidth=2,
                                       #linestyle='-',
                                       alpha=0.6, 
                                       transform=proj_latlon))
    
    # add fire outline
    plt.sca(ax_map)
    plotting.map_fire(ff.data, lats, lons)
    
    
    ax_fp = plt.subplot(3,1,2) # firepower axis
    color_fp = "orange"
    ax_T = ax_fp.twinx() # Temperature
    color_T = "red"
    ax_RH = ax_fp.twinx() # Rel Humid
    color_RH = "magenta"
    ax_ws = plt.subplot(3,1,3, sharex=ax_fp) # Wind speed
    color_ws = 'black'
    color_qc = "darkblue"
    color_wd = "brown"
    ax_ws500 = ax_ws.twinx() # wind speed at 500m
    color_ws500 = 'chocolate'
    
    # offsets for extra y axes, and make just the one spine visible
    ax_RH.spines["right"].set_position(("axes", 1.07))
    plotting.make_patch_spines_invisible(ax_RH)
    ax_RH.spines["right"].set_visible(True)
    
    ## Plot firepower, and PFT on one axis
    plt.sca(ax_fp)
    line_pft, = plt.plot_date(ptimes, pft, '--', 
                             color=color_fp, label='PFT', alpha=0.6)
    line_fp, = plt.plot_date(ftimes, firepower, '-',
                            color=color_fp, label='firepower')
    ## plot temperature on right axis
    plt.sca(ax_T)
    line_T, = plt.plot_date(ctimes, T-273.15, '-',
                           color=color_T, label="Temperature",)
    ## plot RH on further right axis
    plt.sca(ax_RH)
    line_RH, = plt.plot_date(ctimes, RH, '-',
                             color=color_RH, label="RH",)
    plt.plot_date(ctimes, qc_frac, 'o',
                  color=color_qc, 
                  fillstyle='none', 
                  mec=color_qc,
                  mew=1,)
    line_qc, = plt.plot_date(ctimes[qc_heavy], qc_frac[qc_heavy], 'o',
                             color=color_qc,
                             label="Clouds",
                             fillstyle='full')
    plt.plot_date(ctimes[qc_mid], qc_frac[qc_mid], 'o',
                  color=color_qc,
                  fillstyle='bottom')
    
    ## Wind speed, wind speed stdev, wind dir quiver
    plt.sca(ax_ws)
    line_ws, = plt.plot_date(ctimes, ws10, 'o',
                             color=color_ws, label="10m wind speed (+IQR)",)
    line_wsmax, = plt.plot_date(ctimes, ws10_max, '^', 
                                color=color_ws, label="10m max wind speed")
    
    # add quartiles
    plt.fill_between(ctimes, ws10_q3, ws10_q1, alpha=0.6, color='grey')
    
    # normalize windspeed for unit length quivers
    wdnorm = np.sqrt(u10_mean**2 + v10_mean**2)
    # dont show quiver at every single point
    qskip = len(wdnorm)//24 # just show 24 arrows
    plt.quiver(ctimes[::qskip], ws10_q1[::qskip], 
               u10_mean[::qskip]/wdnorm[::qskip], 
               v10_mean[::qskip]/wdnorm[::qskip], 
               pivot='mid', 
               alpha=.9, 
               color=color_ws,
               headwidth=3, headlength=2, headaxislength=2,
               width=.004)
    ## Same again but for 500m wind
    plt.sca(ax_ws500)
    line_ws500, = plt.plot_date(ctimes, ws500, '.',
                                color=color_ws500, 
                                label="500m wind speed (+IQR)",)
    # add quartiles
    plt.fill_between(ctimes, ws500_q3, ws500_q1, alpha=0.4, color=color_ws500)
    
    # normalize windspeed for unit length quivers
    wdnorm500 = np.sqrt(u500_mean**2 + v500_mean**2)
    
    plt.quiver(ctimes[::qskip], ws500_q3[::qskip], 
               u500_mean[::qskip]/wdnorm500[::qskip], 
               v500_mean[::qskip]/wdnorm500[::qskip], 
               pivot='mid', 
               alpha=.7, 
               color=color_wd,
               headwidth=3, headlength=2, headaxislength=2,
               width=.003)
    
    # top row:
    ax_fp.set_ylabel("GWatts")
    ax_fp.yaxis.label.set_color(color_fp)
    ax_T.set_ylabel("Temperature (C)")
    ax_T.yaxis.label.set_color(color_T)
    ax_RH.set_ylabel("RH and rain (frac)")
    ax_RH.yaxis.label.set_color(color_RH)
    
    # bottom row:
    ax_ws.set_ylabel("wind speed (m/s)")
    ax_ws.yaxis.label.set_color(color_ws)
    ax_ws.set_xlabel("Time")
    ax_ws500.set_ylabel("500m wind speed (m/s)")
    ax_ws500.yaxis.label.set_color(color_ws500)
    
    ## Plot periphery
    lines = [line_pft, line_fp, line_T, line_RH, line_qc]
    labels = [l.get_label() for l in lines]
    ax_fp.legend(lines, labels, loc='upper left')
    #plt.suptitle(model_run,y=0.9)
    
    lines2 = [line_ws, line_wsmax, line_ws500]
    labels2 = [l.get_label() for l in lines2]
    ax_ws.legend(lines2,labels2, loc='best')
    
    # format x-ticks date times
    ax_fp.xaxis.set_major_locator(mdates.HourLocator(interval=3))
    ax_fp.xaxis.set_minor_locator(mdates.HourLocator())
    ax_fp.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    
    ax_map.set_title("Surface area-average over time")
    fig.tight_layout(rect=[0,0,.99,1]) # left, bottom, right, top
    fio.save_fig(model_run, _sn_, 'timeseries.png', plt=plt)

if __name__=='__main__':
    
    ## run for all of waroona_run2 datetimes
    #weather_summary_model(model_version='waroona_run3')
    
    ## Run timeseries
    weather_series('waroona_run3')
    
    ## run zoomed in
    #zoom_in = plotting._extents_['sirivans']
    #zoom_in = None # or not
    #weather_summary_model('waroona_run3',zoom_in=zoom_in,HSkip=None, 
    #                      fdtimes=[datetime(2016,1,6,13)])

    print("INFO: weather_summary.py done")

