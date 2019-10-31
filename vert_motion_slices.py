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
from datetime import datetime
import warnings

from utilities import utils,fio,plotting, constants

###
## GLOBALS
###
_sn_ = 'vert_motion_slices'

def vert_motion_slices(qc,w,lh,lat,lon,
                       ff=None,
                       extentname=None,
                       ext='.png',
                       cloud_threshold=constants.cloud_threshold,
                       dpi=300):
    '''
    44i showing vertical motion contourf plots, at different model levels
    Trying to see pyroCB cloud
    '''
    
    # colour bar stuff
    cmap=plotting._cmaps_['verticalvelocity']
    norm=col.SymLogNorm(0.25)
    cbarform=tick.ScalarFormatter()
    contours=np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0]))
    
    # get plot extent
    extent = [lon[0],lon[-1],lat[0],lat[-1]]
    if lon[0] > 130:
        top_mlvl = w.shape[0]-5 # higher for sirivan (where i already cut off the top)
    else:
        top_mlvl = w.shape[0]-20
    plt.close()
    f=plt.figure(figsize=[11,10])
    m_lvls = np.round(np.linspace(0,top_mlvl,16)).astype(int)
    for i in range(16):
        plt.subplot(4,4,i+1)
        
        # top panel is wind speed surface values
        cs,cb=plotting.map_contourf(extent,w[m_lvls[i]],lat,lon,
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
        if extentname is not None:
            plotting.map_add_locations_extent(extentname,hide_text=True)
            
        # Add fire front to first column
        if ff is not None:
            # Add fire outline
            with warnings.catch_warnings():
                # ignore warning when there is no fire yet:
                warnings.simplefilter('ignore')
                plt.contour(lon,lat,np.transpose(ff),np.array([0]), 
                            colors='red',linewidths=1, alpha=0.5 + 0.5*(i==0))
    
    
    # add colour bar
    axes=[0.87, 0.20, 0.04, 0.6] 
    #f.subplots_adjust(top=0.95)
    f.subplots_adjust(wspace=0.001)
    f.subplots_adjust(right=axes[0]-0.03)
    cbar_ax = f.add_axes(axes)
    f.colorbar(cs, cax=cbar_ax, format=cbarform)

def vert_motion_hour(dtime=datetime(2016,1,6,7), model_run='waroona_run1'):
    '''
    create vert motion slices of an hours output from the old run in mika's folder
    '''
    dpi=200
    extentname=model_run.split('_')[0]
    extent=plotting._extents_[extentname]
    ext='.png'
    
    # Read vert motion, clouds
    cubes = fio.read_model_run(model_run, fdtime=dtime, extent=extent)
    w,  = cubes.extract('upward_air_velocity')
    qc, = cubes.extract('qc')
    
    ## fire front
    ff_dtimes = utils.dates_from_iris(w)
    ff, = fio.read_fire(model_run=model_run, dtimes=ff_dtimes, extent=extent, firefront=True)
    
    lh  = w.coord('level_height').points
    lat = w.coord('latitude').points
    lon = w.coord('longitude').points
    
    for i in range(len(ff_dtimes)):
        
        vert_motion_slices(qc[i].data, w[i].data,lh,lat,lon,
                           ff=ff[i].data,
                           extentname=extentname,
                           dpi=dpi, ext=ext)                
        
        stitle = ff_dtimes[i].strftime("%Y %b %d %H:%M (UTC)")
        # Save figure into animation folder with numeric identifier
        plt.suptitle(stitle)
        fio.save_fig(model_run,_sn_,ff_dtimes[i],plt,
                     ext=ext,dpi=dpi)

### RUN THE CODE:
if __name__ == '__main__':
    # set font sizes etc.
    plotting.init_plots()
    
    testing=False
    
    for mr in ['waroona_run1','waroona_old']: #['sirivan_run1','waroona_old', 'waroona_run1']:
        hours = fio.model_outputs[mr]['filedates']
        if testing:
            hours = hours[:3] 
        for dtime in hours:
            vert_motion_hour(dtime, model_run=mr)
    