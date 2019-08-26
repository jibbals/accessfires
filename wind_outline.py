#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 20:45:05 2019
    Winds outline script
@author: jesse
"""

import matplotlib
# don't plot on screen, send straight to file
# this is for lack of NCI display
matplotlib.use('Agg')

# plotting stuff
import matplotlib.colors as col
import matplotlib.pyplot as plt
#import matplotlib.ticker as tick
#import matplotlib.patches as mpatches
import numpy as np
from datetime import datetime,timedelta
import iris # file reading and constraints etc

# local modules
from utilities import plotting, utils, fio


def winds_2panel(s,u,v,w,
                 qc, ff, topog,
                 z,lat,lon,
                 dtime,
                 extentname='waroona',
                 transect=1, 
                 vectorskip=9,
                 quiverscale=60,
                 ext='.png'
                ):
    '''
    311 Plot showing contour map and wind speed, along with near sites and transect
    312 plot showing vertical motion along transect
    313 plot showing wind speed along transect
    INPUTS:
        dictionary with topog, vert motion, wind speed, wind direction, zth,lat,lon, time
        timestep: data will have a time dimension
        extentname = { 'waroona' | 'sirivan' } choice of two fire locations
        transect = int from 1 to 6 for transect choice
        vectorskip reduces quiver density (may need to fiddle)
        quiverscale changes how long the arrows are (also may need to fiddle)
        ext is the plot extension { '.png' | '.eps' }
    '''
    
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
    plt.subplot(3,1,1)
    
    # top panel is topography
    plotting.map_topography(extent,topog,lat,lon)
    plt.title('Topography, winds')
    
    # Add fire front contour
    if ff is not None:
        plt.contour(lon,lat,np.transpose(ff),np.array([0]), 
                    colors='red')
    
    # start to end x=[lon0,lon1], y=[lat0, lat1]
    plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
             linewidth=2, 
             marker='X', markersize=7,markerfacecolor='white')
    
    # add nearby towns
    textcolor='k'
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
    plt.get_cmap(plotting._cmaps_['windspeed'])
    plt.quiver(lon[skip[1]],lat[skip[0]],u[0][skip],v[0][skip], scale=quiverscale)
               #color=cmap(norm(s[skip])), 
    
    
    ## Second row is transect plots
    plt.subplot(3,1,2)
    plotting.transect_w(w,z, lat, lon,start,end,topog=topog)
    plt.ylabel('height (m)')
    #plt.xlabel('transect')
    
    ax3 = plt.subplot(3,1,3)
    trs, trx, trz = plotting.transect_s(s,z,lat,lon,start,end,topog=topog)
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
    
    # Save figure into animation folder with numeric identifier
    plt.suptitle(stitle)
    print("INFO: Saving figure:",pname)
    plt.savefig(pname)
    plt.close()
    return pname


def waroona_wind_loop(dtime):
    '''
    Loop over transects over waroona and make the wind outline figures
    '''
    um_hour=datetime(dtime.year,dtime.month,dtime.day,dtime.hour)
    if um_hour < datetime(2016,1,6,15):
        ffpath = 'data/waroona_fire/firefront.CSIRO_24h.20160105T1500Z.nc'
    elif um_hour < datetime(2016,1,7,15):
        ffpath = 'data/waroona_fire/firefront.CSIRO_24h.20160106T1500Z.nc'
    else:
        ffpath = 'data/waroona_fire/firefront.CSIRO_24h.20160107T1500Z.nc'
    
    extentname='waroona'
    
    # Constraints on dimensions (save ram and reading time)
    West,East,South,North = plotting._extents_['waroona']
    constr_z = iris.Constraint(model_level_number=lambda cell: cell < 90)
    constr_lons = iris.Constraint(longitude = lambda cell: West <= cell <= East )
    constr_lats = iris.Constraint(latitude = lambda cell: South <= cell <= North )
    
    ## Model level heights and topography don't depend on time
    # metres
    zth, = iris.load('data/waroona/umnsaa_2016010515_mdl_th1.nc',
                     ['height_above_reference_ellipsoid' &
                      constr_z & 
                      constr_lats & 
                      constr_lons])
    
    topog, = fio.read_nc_iris('data/waroona/umnsaa_2016010515_slv.nc',
                              constraints = 'surface_altitude'  & 
                              constr_lats & 
                              constr_lons)
    
    
    # Read the cubes
    slv,ro1,th1,th2 = fio.read_waroona_iris(dtime=um_hour, 
                                            constraints = [constr_z &
                                                           constr_lons &
                                                           constr_lats])
    
    # wind speeds need to be interpolated onto non-staggered latlons
    p, u1, v1 = ro1.extract(['air_pressure','x_wind','y_wind'])
    # DESTAGGER u and v using iris interpolate
    # u1: [t,z, lat, lon1]
    # v1: [t,z, lat1, lon]  # put these both onto [t,z,lat,lon]
    u = u1.interpolate([('longitude',p.coord('longitude').points)],
                       iris.analysis.Linear())
    v = v1.interpolate([('latitude',p.coord('latitude').points)],
                       iris.analysis.Linear())
    
    
    # Get wind speed cube using hypotenuse of u,v (I think this is the first action that actually reads any file data)
    s = iris.analysis.maths.apply_ufunc(np.hypot,u,v)

    # Cloud parameters and upward wind velocity    
    qc1,qc2 = th2.extract(['mass_fraction_of_cloud_ice_in_air','mass_fraction_of_cloud_liquid_water_in_air'])    
    qc = ( qc1+qc2 )*1000 # change to g/kg    
    w, = th1.extract('upward_air_velocity')
    
    ## fire front
    # read 6 time steps:
    ff_dtimes = np.array([um_hour + timedelta(hours=x/60.) for x in range(10,61,10)])
    ff = fio.read_fire_front(ffpath,dtimes=ff_dtimes)
    ff = ff.extract(constr_lats & constr_lons) # subset lats,lons
    
    # datetime of outputs
    tdim = p.coord('time')
    d0 = datetime.strptime(str(tdim.units),'hours since %Y-%m-%d %H:%M:00')
    timesteps = utils.date_from_gregorian(tdim.points, d0=d0)
    
    lat,lon = w.coord('latitude').points, w.coord('longitude').points
    
    # also loop over different transects
    for i in range(5):
        for tstep in range(len(timesteps)):
            print("DEBUG: ", s.shape, qc.shape, ff.shape, topog.shape, zth.shape, timesteps[tstep])
            winds_2panel(s[tstep].data, u[tstep].data, v[tstep].data, w[tstep].data,
                         qc[tstep].data, ff[tstep].data, topog.data,
                         zth.data, lat,lon,
                         dtime=timesteps[tstep],
                         extentname=extentname,
                         transect=i+1)
    
    

#########################################################################
#########################################################################
#########################################################################
    
    
if __name__ == '__main__':
    
    print("INFO: testing wind_outline.py")
    waroona_wind_loop(datetime(2016,1,5,15))
    #for dtime in [ datetime(2016,1,6,7) + timedelta(hours=x) for x in range(4) ]:
        
    #    waroona_wind_loop(dtime)