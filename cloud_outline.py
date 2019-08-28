# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 15:58:21 2019
    Create cloud summary plots
@author: jgreensl
"""

import matplotlib
# don't plot on screen, send straight to file
# this is for lack of NCI display
matplotlib.use('Agg',warn=False)

# plotting stuff
#import matplotlib.colors as col
import matplotlib.pyplot as plt
#import matplotlib.ticker as tick
#import matplotlib.patches as mpatches
import numpy as np
from datetime import datetime,timedelta
# ignore some warnings
import warnings

# local modules
from utilities import plotting, utils, fio, constants


def clouds_2panel(topog,s,u,v,
                  qc,w,
                  ff,
                  z,lat,lon,
                  dtime,
                  extentname='waroona',
                  transect=1, 
                  vectorskip=13,
                  quiverscale=60,
                  ztop=10000,
                  ext='.png',
                  cloud_threshold=constants.cloud_threshold,
                  dpi=400,
                  ):
    '''
    311 Plot showing windspeed map and wind speed, along with near sites and transect
    312 relative humidity along transect
    313 water and ice content along transect
    INPUTS:
        cubes need to be passed in
        extentname = { 'waroona' | 'sirivan' } choice of two fire locations
        transect = int from 1 to 6 for transect choice (see figures/transects.png)
        vectorskip reduces quiver density (may need to fiddle)
        quiverscale changes how long the arrows are (also may need to fiddle)
        dtime is datetime of output 
        ext is the plot extension { '.png' | '.eps' }
    '''
    
    ## Plot data, inputs will be [[z],lat,lon]
    
    # datetime timestamp for file,title
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
    plt.subplot(3,1,1)
    
    # top panel is wind speed surface values
    plotting.map_contourf(extent,s,lat,lon,
                          clevs = np.linspace(0,15,31),
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
    plt.quiver(lon[skip[1]],lat[skip[0]],u[skip],v[skip], scale=quiverscale)
    
    # Add fire outline
    plt.contour(lon,lat,np.transpose(ff),np.array([0]), 
                colors='red')
    
    ## Second row is transect plots
    plt.subplot(3,1,2)
    #plotting.transect(theta,z,lat,lon,start,end,topog=topog,
    #                  cmap='plasma',
    #                  contours=np.linspace(290,400,111),
    #                  ztop=ztop)
    wslice, xslice, zslice = plotting.transect_w(w,z,lat,lon,start,end,topog=topog,
                                                 npoints=100,lines=None,
                                                 ztop=ztop)
    plt.title("Vertical motion (m/s)")
    
    qcslice = utils.cross_section(qc,lat,lon,start,end, npoints=100)
    with warnings.catch_warnings():
        # ignore warning when there are no clouds:
        warnings.simplefilter('ignore')
        plt.contour(xslice,zslice,qcslice,np.array([cloud_threshold]),colors='teal')
    
    plt.ylabel('height (m)')
    plt.xlabel('')
    
    plt.subplot(3,1,3)
    plotting.transect_qc(qc,z,lat,lon,start,end,topog=topog,
                        ztop=ztop,)
    # Show transect start and end
    xticks,xlabels = plotting.transect_ticks_labels(start,end)
    plt.xticks(xticks[0::2],xlabels[0::2]) # show transect start and end
    
    # Save figure into animation folder with numeric identifier
    plt.suptitle(stitle)
    print("INFO: Saving figure:",pname)
    plt.savefig(pname,dpi=dpi)
    plt.close()
    return pname


def waroona_cloud_loop(dtime):
    '''
    make an hours worth of clouds_2panel plots starting at argument dtime
    First get iris cubes from each of the data files we will read,
        subset the data before reading it to save ram and run faster
        also read fire front at matching times
    then send all the data to plotting method
    '''
    um_hour=datetime(dtime.year,dtime.month,dtime.day,dtime.hour)
    extentname='waroona'
    extent=plotting._extents_[extentname]
    um_hour=datetime(dtime.year,dtime.month,dtime.day,dtime.hour)
    
    # Read the cubes
    slv,ro1,th1,th2 = fio.read_waroona(dtime=um_hour, extent=extent,
                                            add_winds=True, add_theta=True)
    
    zro, = ro1.extract('z_ro')
    
    topog, = slv.extract('topog')
    
    u,v,s = ro1.extract(['u','v','s'])

    qc, = th2.extract('qc')
    
    w, = th1.extract('upward_air_velocity')
    
    ## fire front
    # read 6 time steps:
    ff_dtimes = np.array([um_hour + timedelta(hours=x/60.) for x in range(10,61,10)])
    ff, = fio.read_fire(dtimes=ff_dtimes,extent=extent, firefront=True)
    
    sh,Ta = th1.extract(['specific_humidity','air_temperature'])
    
    theta, = th1.extract('potential_temperature')
    
    lat,lon = w.coord('latitude').points, w.coord('longitude').points
    
    # datetime of outputs
    tdim = w.coord('time')
    d0 = datetime.strptime(str(tdim.units),'hours since %Y-%m-%d %H:%M:00')
    timesteps = utils.date_from_gregorian(tdim.points, d0=d0)
    
    # also loop over different transects
    for i_transect in np.arange(1,6.5,1, dtype=int):
        for tstep in range(len(timesteps)):
            clouds_2panel(topog.data, s[tstep,0].data, u[tstep,0].data, v[tstep,0].data,
                          qc[tstep].data, w[tstep].data, 
                          ff[tstep].data,
                          zro.data, lat, lon,
                          dtime=timesteps[tstep],
                          extentname=extentname,
                          transect=i_transect)
            
if __name__ == '__main__':
    
    print("INFO: testing cloud_outline.py")
    #waroona_cloud_loop(datetime(2016,1,5,15))
    for dtime in [ datetime(2016,1,6,7) + timedelta(hours=x) for x in range(3) ]:
        
        waroona_cloud_loop(dtime)

