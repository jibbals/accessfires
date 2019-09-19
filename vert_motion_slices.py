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

###
## GLOBALS
###
_sn_ = 'vert_motion_slices'

def vert_motion_slices(qc,w,lh,lat,lon,dtime, 
                       ff=None,
                       model_run='waroona_run1',
                       ext='.png',
                       cloud_threshold=constants.cloud_threshold,
                       dpi=300):
    '''
    44i showing vertical motion contourf plots, at different model levels
    Trying to see pyroCB cloud
    '''
    
    extentname=model_run.split('_')[0]
    # datetime timestamp for file,title
    stitle = dtime.strftime("%Y %b %d %H:%M (UTC)")
    
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
    
    fio.save_fig(model_run,_sn_,dtime,plt,
                 ext=ext,dpi=dpi)
    

def vert_motion_hour(dtime=datetime(2016,1,6,7), model_run='waroona_run1'):
    '''
    create vert motion slices of an hours output from the old run in mika's folder
    '''
    dpi=200
    extentname=model_run.split('_')[0]
    extent=plotting._extents_[extentname]
    ext='.png'
    
    # Read vert motion, clouds
    if model_run=='waroona_old':
        cubes = fio.read_waroona_pcfile(dtime,extent=extent)
        w,  = cubes.extract('upward_air_velocity')
        qc, = cubes.extract('qc')
    elif model_run=='waroona_run1':
        _,_,th1cubes,th2cubes = fio.read_waroona(dtime,extent=extent)
        w,  = th1cubes.extract('upward_air_velocity')
        qc, = th2cubes.extract('qc')
    else:
        assert False, '%s not yet implemented for run %s'%(_sn_,model_run)
    
    ## fire front
    ff_dtimes = utils.dates_from_iris(w)
    ff1, = fio.read_fire(dtimes=ff_dtimes, extent=extent, firefront=True)
    ff=ff1
    if model_run=='waroona_old':
        ff=ff1.interpolate([('longitude',w.coord('longitude').points), 
                            ('latitude',w.coord('latitude').points)],
                           iris.analysis.Linear())
    
    lh  = w.coord('level_height').points
    lat = w.coord('latitude').points
    lon = w.coord('longitude').points
    
    for i in range(len(ff_dtimes)):
        subtime = ff_dtimes[i]
        
        vert_motion_slices(qc[i].data, w[i].data,lh,lat,lon,
                           ff=ff[i].data,
                           dtime=subtime,
                           model_run=model_run,
                           dpi=dpi, ext=ext)                

### RUN THE CODE:
if __name__ == '__main__':
    ### pyrocb utc time window:
    run1_hours = [datetime(2016,1,6,7) + timedelta(hours=x) for x in range(2)]
    old_hours = [datetime(2016,1,6,4) + timedelta(hours=x) for x in range(5)]
    
    for dtime in old_hours:
        vert_motion_hour(dtime, model_run='waroona_old')
    for dtime in run1_hours:
        vert_motion_hour(dtime, model_run='waroona_run1')
    