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

from utilities import fio,plotting, constants


def vert_motion_slices(qc,w,lh,lat,lon,dtime, 
                       extentname='waroona',
                       ext='.png',
                       cloud_threshold=constants.cloud_threshold,
                       dpi=400):
    '''
    43i showing vertical motion contourf plots, at different model levels
    Trying to see pyroCB cloud
    '''
    
    # datetime timestamp for file,title
    dstamp = dtime.strftime("%Y%m%d%H%M")
    stitle = dtime.strftime("%Y %b %d %H:%M (UTC)")
    
    # figure name and location
    pname="figures/%s/vert_motion_slices/fig_%s%s"%(extentname,dstamp,ext)
    
    
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
    f=plt.figure(figsize=[10,10])
    # model levels to plot
    m_lvls = [30,60,72,  20,50,69,  10,40,66,  5,35,63] 
    for i in range(12):
        plt.subplot(4,3,i+1)
        
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



### RUN THE CODE:
if __name__ == '__main__':
    dpi=400
    extentname='waroona'
    extent=plotting._extents_[extentname]
    ext='.png'
    
    ### pyrocb utc time window:
    pyrocb_hours = [datetime(2016,1,6,7) + timedelta(hours=x) for x in range(6)]
    #test hours
    #pyrocb_hours = [datetime(2016,1,5,15)]
    
    for dtime in pyrocb_hours:
        
        # Read vert motion
        _,_,th1cubes,th2cubes = fio.read_waroona(dtime,extent=extent)
        w,  = th1cubes.extract('upward_air_velocity')
        qc, = th2cubes.extract('qc')
        lh  = w.coord('level_height').points
        lat = w.coord('latitude').points
        lon = w.coord('longitude').points
        
        for i in range(6):
            subtime = dtime+timedelta(seconds=600*(i+1))
            vert_motion_slices(qc[i].data, w[i].data,lh,lat,lon,dtime=subtime)
