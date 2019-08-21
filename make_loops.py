
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  20 2019

  Script to make the plots shown in animators.ipynb
  
@author: jesse
"""

import matplotlib
matplotlib.use('Agg')# don't plot on screen, send straight to file
# this is for NCI display issue
import matplotlib.pyplot as plt
plt.ioff() # fix problem with ploting datetimes (don't know why this works)

import matplotlib as mpl
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
from datetime import datetime,timedelta

from multiprocessing import Pool
from time import sleep

from utilities import plotting, utils, fio


def winds_2panel(data,tstep, 
                 extentname='waroona',
                 transect=1, 
                 vectorskip=14,
                 quiverscale=60,
                 dtime=None,
                 ext='.png'
                ):
    '''
    211 Plot showing contour map and wind speed, along with near sites and transect
    223 plot showing vertical motion along transect
    224 plot showing wind speed along transect
    INPUTS:
        topography, vert motion, wind speed, z,lat,lon, 
        extentname = { 'waroona' | 'sirivan' } choice of two fire locations
        transect = int from 1 to 6 for transect choice
        vectorskip reduces quiver density (may need to fiddle)
        quiverscale changes how long the arrows are (also may need to fiddle)
        dtime is datetime of output 
        ext is the plot extension { '.png' | '.eps' }
    '''
    topog=data['topog']
    w=data['upward_air_velocity'][tstep]
    s=data['wind_speed'][tstep]
    u=data['x_wind_destaggered'][tstep]
    v=data['y_wind_destaggered'][tstep]
    z=data['zth'][tstep]
    lat=data['latitude']
    lon=data['longitude']
    
    
    # datetime timestamp for file,title
    if dtime is None:
        dstamp = "20YYMMDDHHMM"
        stitle = dstamp
    else:
        dstamp = dtime.strftime("%Y%m%d%H%M")
        stitle = dtime.strftime("%Y %b %d %H:%M (UTC)")
    # figure name and location
    pname="figures/%s/winds_outline_X%d/fig_%s%s"%(extentname,transect,dstamp,ext)
    
    # set font sizes
    plotting.init_plots()
    
    # get plot extent, and transect
    extent = plotting._extents_[extentname]
    start,end = plotting._transects_["%s%d"%(extentname,transect)]
    
    plt.figure(figsize=[7,10])
    ax1 = plt.subplot(3,1,1)
    
    # top panel is topography
    plotting.map_topography(extent,topog,lat,lon)
    plt.title('Topography, winds')
    
    # start to end x=[lon0,lon1], y=[lat0, lat1]
    plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
             linewidth=2, 
             marker='X', markersize=7,markerfacecolor='white')
    
    # add nearby towns
    textcolor='lightgrey'
    if extentname == 'waroona':
        plotting.map_add_locations(['waroona','yarloop'], 
                                   text=['Waroona', 'Yarloop'], 
                                   textcolor=textcolor)
        # add fire ignition
        plotting.map_add_locations(['fire_waroona'],
                                   text = ['Fire ignition'], 
                                   color='r', marker='*', 
                                   textcolor=textcolor)
        # add pyroCB
    else:
        plotting.map_add_locations(['sirivan','uarbry'], 
                                   text=['Sir Ivan','Uarbry'],
                                   dx=[-.02,.05], dy =[-.015,-.03],
                                   textcolor=textcolor)
        # add fire ignition
        plotting.map_add_locations(['fire_sirivan'],
                                   text = ['Fire ignition'], dx=.05,
                                   color='r', marker='*', 
                                   textcolor=textcolor)
        # add pyroCB

    
    # Add vectors for winds
    # just surface, and one every 10 points to reduce density
    skip = (slice(None,None,vectorskip),slice(None,None,vectorskip))
    ## colour the arrows
    # map wind speed to colour map domain [0, 1]
    norm = col.Normalize()
    norm.autoscale(s[skip])
    cmap = plt.get_cmap(plotting._cmaps_['windspeed'])
    plt.quiver(lon[skip[1]],lat[skip[0]],u[0][skip],v[0][skip], scale=quiverscale)
               #color=cmap(norm(s[skip])), 
    
    
    ## Second row is transect plots
    ax2 = plt.subplot(3,1,2)
    plotting.transect_w(w,z, lat, lon,start,end,topog=topog)
    plt.ylabel('height (m)')
    #plt.xlabel('transect')
    
    ax3 = plt.subplot(3,1,3)
    contours=np.arange(0,37,3)
    trs, trx, trz = plotting.transect_s(s,z,lat,lon,start,end,topog=topog,
                                        contours=contours,lines=contours)
    #plt.yticks([])
    xticks,xlabels = plotting.transect_ticks_labels(start,end)
    plt.xticks(xticks[0::2],xlabels[0::2]) # show transect start and end
    
    # Annotate max wind speed
    # only care about lower troposphere
    # 60 levels is about 2800m, 80 levels is about 4700m, 70 levels : 3700m
    upto = 70 
    mloc = np.unravel_index(np.argmax(trs[:upto,:],axis=None),trs[:upto,:].shape)
    note="max = %5.1f m/s\n  (at %4.0f m) "%(trs[:upto,:][mloc], trz[:upto,:][mloc])
    trans = ax3.get_xaxis_transform() # x in data untis, y in axes fraction
    ann = ax3.annotate(note, fontsize=15,
                       xy=(0.33, -0.2 ), xycoords=trans)
    
    # Save figure into animation folder with numeric identifier
    plt.suptitle(stitle)
    print("INFO: Saving figure:",pname)
    plt.savefig(pname)
    plt.close()
    return pname


def clouds_2panel(data, tstep,
                 extentname='waroona',
                 transect=1, 
                 vectorskip=14,
                 quiverscale=60,
                 ztop=13000,
                 dtime=None,
                 ext='.png'
                ):
    '''
    211 Plot showing windspeed map and wind speed, along with near sites and transect
    223 vertical motion along transect
    224 water and ice content along transect
    INPUTS:
        topography, vert motion, wind speed, u,v, qc, z,lat,lon, 
        extentname = { 'waroona' | 'sirivan' } choice of two fire locations
        transect = int from 1 to 6 for transect choice
        vectorskip reduces quiver density (may need to fiddle)
        quiverscale changes how long the arrows are (also may need to fiddle)
        dtime is datetime of output 
        ext is the plot extension { '.png' | '.eps' }
    '''
    
    topog=data['topog']
    #w=data['upward_air_velocity'][tstep]
    sh=data['specific_humidity'][tstep]
    s=data['wind_speed'][tstep]
    u=data['x_wind_destaggered'][tstep]
    v=data['y_wind_destaggered'][tstep]
    z=data['zth'][tstep]
    lat=data['latitude']
    lon=data['longitude']
    qc=data['qc'][tstep]
    
    # datetime timestamp for file,title
    if dtime is None:
        dstamp = "YYYYMMDDHHMM"
        stitle = dstamp
    else:
        dstamp = dtime.strftime("%Y%m%d%H%M")
        stitle = dtime.strftime("%Y %b %d %H:%M (UTC)")
    # figure name and location
    pname="figures/%s/clouds_outline_X%d/fig_%s%s"%(extentname,transect,dstamp,ext)
    
    # set font sizes
    plotting.init_plots()
    
    # get plot extent, and transect
    extent = plotting._extents_[extentname]
    start,end = plotting._transects_["%s%d"%(extentname,transect)]
    
    plt.figure(figsize=[7,10])
    ax1 = plt.subplot(3,1,1)
    
    # top panel is wind speed surface values
    plotting.map_contourf(extent,s[0],lat,lon,
                          cmap=plotting._cmaps_['windspeed'],
                          clabel='m/s')
    plt.title('Horizontal wind speed')
    

    # start to end x=[lon0,lon1], y=[lat0, lat1]
    plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
             linewidth=2, 
             marker='X', markersize=7,markerfacecolor='white')
    
    # add nearby towns
    if extentname == 'waroona':
        plotting.map_add_locations(['waroona','yarloop'], 
                                   text=['Waroona', 'Yarloop'], 
                                   textcolor='k')
        # add fire ignition
        plotting.map_add_locations(['fire_waroona'],
                                   text = ['Fire ignition'], 
                                   color='r', marker='*', 
                                   textcolor='k')
        # add pyroCB
    else:
        plotting.map_add_locations(['sirivan','uarbry'], 
                                   text=['Sir Ivan','Uarbry'],
                                   dx=[-.02,.05], dy =[-.015,-.03],
                                   textcolor='k')
        # add fire ignition
        plotting.map_add_locations(['fire_sirivan'],
                                   text = ['Fire ignition'], dx=.05,
                                   color='r', marker='*', 
                                   textcolor='k')
        # add pyroCB

    
    # Add vectors for winds
    # just surface, and one every N points to reduce density
    skip = (slice(None,None,vectorskip),slice(None,None,vectorskip))
    #mlon,mlat = np.meshgrid(lon,lat)
    plt.quiver(lon[skip[1]],lat[skip[0]],u[0][skip],v[0][skip], scale=quiverscale)
    
    
    ## Second row is transect plots
    ax2 = plt.subplot(3,1,2)
    plotting.transect(sh,z,lat,lon,start,end,topog=topog,
                      cmap='plasma',
                      ztop=ztop)
    plt.title("Specific humidity")
    #plotting.transect_w(w,z, lat, lon,start,end,topog=topog)
    
    plt.ylabel('height (m)')
    plt.xlabel('')
    
    ax3 = plt.subplot(3,1,3)
    plotting.transect_qc(qc,z,lat,lon,start,end,topog=topog,
                        ztop=ztop,)
    # Show transect start and end
    xticks,xlabels = plotting.transect_ticks_labels(start,end)
    plt.xticks(xticks[0::2],xlabels[0::2]) # show transect start and end
    
    # Save figure into animation folder with numeric identifier
    plt.suptitle(stitle)
    print("INFO: Saving figure:",pname)
    plt.savefig(pname)
    plt.close()
    return pname

def waroona_cloud_loop(dtime):
    '''
    make an hours worth of clouds_2panel plots
    '''
    # Waroona loop
    extentname='waroona'
    vectorskip=9

    # List of hours for which we have output.. 0515 - 0614
    #date_list = [datetime(2016,1,5,15) + timedelta(hours=x) for x in range(24)]
    # first hour is already done, subsequent hours have some missing fields!
    topography,latt,lont = fio.read_topog('data/waroona/topog.nc')

    # Read the files for this hour
    waroona_outputs = fio.read_waroona(dtime)
    slv, ro1, th1, th2 = waroona_outputs
    
    # grab the outputs desired
    waroona={}
    waroona['wind_speed'] = ro1['wind_speed']
    waroona['x_wind_destaggered'] = ro1['x_wind_destaggered']
    waroona['y_wind_destaggered'] = ro1['y_wind_destaggered']
    waroona['topog'] = topography
    waroona['zth'] = th1['zth']
    waroona['latitude'] = slv['latitude']
    waroona['longitude'] = slv['longitude']
    waroona['specific_humidity'] = th1['specific_humidity']
    waroona['qc'] = th2['mass_fraction_of_cloud_ice_in_air'] + \
                    th2['mass_fraction_of_cloud_liquid_water_in_air']
    #waroona['upward_air_velocity'] = th1['upward_air_velocity']
    
    
    # datetime of outputs
    timesteps = utils.date_from_gregorian(th1['time'])
    
    # also loop over different transects
    for i_transect in np.arange(1,6.5,1, dtype=int):
        for tstep in range(len(timesteps)):
            pname = clouds_2panel(waroona,tstep=tstep,
                                 dtime=timesteps[tstep],
                                 extentname=extentname,
                                 vectorskip=vectorskip, 
                                 transect=i_transect)
    
    # Save ram hopefully
    del waroona, slv, ro1, th1, th2, waroona_outputs

def waroona_wind_loop(dtime):
    '''
    Loop over transects over waroona and make the wind outline figures
    '''
    # Waroona loop
    extentname='waroona'
    vectorskip=9

    # List of hours for which we have output.. 0515 - 0614
    #date_list = [datetime(2016,1,5,15) + timedelta(hours=x) for x in range(24)]
    # first hour is already done, subsequent hours have some missing fields!
    topography,latt,lont = fio.read_topog('data/waroona/topog.nc')

    # Read the files for this hour
    waroona_outputs = fio.read_waroona(dtime,th2=False) # don't need cloud stuff for this plot
    slv, ro1, th1, th2 = waroona_outputs
    
    # grab the outputs desired
    waroona={}
    waroona['wind_speed'] = ro1['wind_speed']
    waroona['x_wind_destaggered'] = ro1['x_wind_destaggered']
    waroona['y_wind_destaggered'] = ro1['y_wind_destaggered']
    waroona['topog'] = topography
    waroona['zth'] = th1['zth']
    waroona['latitude'] = slv['latitude']
    waroona['longitude'] = slv['longitude']
    waroona['upward_air_velocity'] = th1['upward_air_velocity']
    
    
    # datetime of outputs
    timesteps = utils.date_from_gregorian(th1['time'])
    
    # also loop over different transects
    for i_transect in np.arange(1,6.5,1, dtype=int):
        for tstep in range(len(timesteps)):
            pname = winds_2panel(waroona,tstep=tstep,
                                 dtime=timesteps[tstep],
                                 extentname=extentname,
                                 vectorskip=vectorskip, 
                                 transect=i_transect)
    
    # Save ram hopefully
    del waroona, slv, ro1, th1, th2, waroona_outputs



if __name__=='__main__':
    
    # TODO make input arguments for each loop being created
    
    # First 24 hours:
    day1_dtimes = [datetime(2016,1,5,15) + timedelta(hours=x) for x in range(24)]
    still_need_times = [datetime(2016,1,6,7) + timedelta(hours=x) for x in range(8)]
    first_10_hrs = day1_dtimes[:10]
    
    for dtime in first_10_hrs:
        waroona_cloud_loop(dtime)
    
    #for dtime in still_need_times:
    #    waroona_wind_loop(dtime)
    
    #with Pool(processes=2) as pool:
        
        ## Send each datetime to the process pool
        #pool.map( waroona_wind_loop, day1_dtimes )
        #pool.map( waroona_cloud_loop, day1_dtimes )
        
