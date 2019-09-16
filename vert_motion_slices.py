# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 15:22:49 2019
    Show vertical motion slices
@author: jgreensl
"""
import matplotlib
matplotlib.use("Agg", warn=False)

import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.ticker as tick

import numpy as np
from datetime import datetime, timedelta
import warnings
import iris

from utilities import utils,fio,plotting, constants


def vert_motion_slices(qc,w,lh,lat,lon,dtime, 
                       ff=None,
                       extentname='waroona',
                       ext='.png',
                       folder='vert_motion_slices',
                       cloud_threshold=constants.cloud_threshold,
                       dpi=400):
    '''
    44i showing vertical motion contourf plots, at different model levels
    Trying to see pyroCB cloud
    '''
    
    # datetime timestamp for file,title
    dstamp = dtime.strftime("%Y%m%d%H%M")
    stitle = dtime.strftime("%Y %b %d %H:%M (UTC)")
    
    # figure name and location
    pname="figures/%s/%s/fig_%s%s"%(extentname,folder,dstamp,ext)
    
    
    # set font sizes etc.
    plotting.init_plots()
    
    # colour bar stuff
    cmap=plotting._cmaps_['verticalvelocity']
    norm=col.SymLogNorm(0.25)
    cbarform=tick.ScalarFormatter()
    contours=np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0]))
    
    # get plot extent
    extent = plotting._extents_[extentname]
    
    plt.close()
    f=plt.figure(figsize=[11,10])
    # model levels to plot
    m_lvls = [30,60,72,90,  20,50,69,85,  10,40,66,80,  5,35,63,76] 
    for i in range(16):
        plt.subplot(4,4,i+1)
        
        # top panel is wind speed surface values
        cs=plotting.map_contourf(extent,w[m_lvls[i]],lat,lon,
                                 clevs = contours,
                                 cmap=cmap, norm=norm,
                                 cbar=False,)#clabel='m/s', cbarform=cbarform)
        
        plt.title("~ %5.0f m"%lh[m_lvls[i]])
        
        ## Add contour where clouds occur
        with warnings.catch_warnings():
            # ignore warning when there are no clouds:
            warnings.simplefilter('ignore')
            plt.contour(lon,lat,qc[m_lvls[i]],np.array([cloud_threshold]), 
                        colors='k')
        
        # add nearby towns
        if extentname == 'waroona':
            plotting.map_add_locations(['waroona','yarloop'], 
                                       text=['', ''], 
                                       textcolor='k')
            # add fire ignition
            plotting.map_add_locations(['fire_waroona'],
                                       text = [''], 
                                       color='r', marker='*', 
                                       textcolor='k')
            # add pyroCB
        # Add fire front to first column
        if ff is not None:
            # Add fire outline
            with warnings.catch_warnings():
                # ignore warning when there is no fire yet:
                warnings.simplefilter('ignore')
                plt.contour(lon,lat,np.transpose(ff),np.array([0]), 
                            colors='red',linewidths=1, alpha=0.5 + 0.5*(i==12))
    
    # Save figure into animation folder with numeric identifier
    plt.suptitle(stitle)
    
    # add colour bar
    axes=[0.87, 0.20, 0.04, 0.6] 
    #f.subplots_adjust(top=0.95)
    f.subplots_adjust(wspace=0.001)
    f.subplots_adjust(right=axes[0]-0.03)
    cbar_ax = f.add_axes(axes)
    f.colorbar(cs, cax=cbar_ax, format=cbarform)
    
    print("INFO: Saving figure:",pname)
    plt.savefig(pname,dpi=dpi)
    plt.close()


def vert_motion_hour(dtime=datetime(2016,1,6,7), old=False):
    '''
    create vert motion slices of an hours output from the old run in mika's folder
    '''
    dpi=400
    extentname='waroona'
    extent=plotting._extents_[extentname]
    ext='.png'
    
    # Read vert motion, clouds
    if old:
        cubes = fio.read_waroona_pcfile(dtime,extent=extent)
        w,  = cubes.extract('upward_air_velocity')
        qc, = cubes.extract('qc')
    else:
        _,_,th1cubes,th2cubes = fio.read_waroona(dtime,extent=extent)
        w,  = th1cubes.extract('upward_air_velocity')
        qc, = th2cubes.extract('qc')
    lh  = w.coord('level_height').points
    lat = w.coord('latitude').points
    lon = w.coord('longitude').points

    ## fire front
    ff_dtimes = utils.dates_from_iris(w)
    ff1, = fio.read_fire(dtimes=ff_dtimes, extent=extent, firefront=True)
    
    # interpolat ff onto old lats and lons    
    ff=ff1
    if old:
        ff=ff1.interpolate([('longitude',w.coord('longitude').points), 
                            ('latitude',w.coord('latitude').points)],
                           iris.analysis.Linear()) 
    for i in range(len(ff_dtimes)):
        subtime = ff_dtimes[i]
        
        vert_motion_slices(qc[i].data, w[i].data,lh,lat,lon,
                           ff=ff[i].data,
                           dtime=subtime,
                           folder='vert_motion_slices%s'%(['','_old'][old]),
                           dpi=dpi, ext=ext)                

### RUN THE CODE:
if __name__ == '__main__':
    ### pyrocb utc time window:
    pyrocb_hours = [datetime(2016,1,6,7) + timedelta(hours=x) for x in range(6)]
    #test hours
    pyrocb_hours = [datetime(2016,1,6,7)]
    
    #for dtime in pyrocb_hours:
    #    plot_hour(dtime,old=True)
    
    #for dtime in [ datetime(2016,1,6,7) + timedelta(hours=x) for x in range(2) ]:
    for dtime in [ datetime(2016,1,5,15) + timedelta(hours=x) for x in range(8) ]:
        vert_motion_hour(dtime,old=True)
    for dtime in [ datetime(2016,1,6,4) + timedelta(hours=x) for x in range(5) ]:
        vert_motion_hour(dtime,old=True)
    for dtime in [ datetime(2016,1,5,15) + timedelta(hours=x) for x in range(24) ]:
        vert_motion_hour(dtime, old=False)
