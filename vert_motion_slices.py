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
from datetime import datetime,timedelta
import warnings

from utilities import utils,fio,plotting, constants

###
## GLOBALS
###
_sn_ = 'vert_motion_slices'

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
    norm=col.SymLogNorm(0.25)
    cbarform=tick.ScalarFormatter()
    contours=np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0]))
    
    # get plot extent
    extent = [lon[0],lon[-1],lat[0],lat[-1]]
    plt.close()
    f=plt.figure(figsize=[11,10])
    m_lvls = np.array([10, 20, 32,  40,  48,  56,  64,  72,  80,  88,  96,
       104, 108, 112, 116, 120],np.int)
    nrows,ncols=4,4
    if three_by_three:
        nrows,ncols=3,3
        m_lvls = m_lvls[np.array([0,3,6,8,10,12,13,14,15])]
    #m_lvls = np.round(np.linspace(0,120,16)).astype(int)
    for i in range(len(m_lvls)):
        plt.subplot(nrows,ncols,i+1)
        
        # top panel is wind speed surface values
        cs,cb=plotting.map_contourf(extent,w[m_lvls[i]],lat,lon,
                                    levels = contours,
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
            
        # Add fire front
        plotting.map_fire(ff,lat,lon,
                          colors='red',
                          linewidths=1, 
                          alpha=0.5 + 0.5*(i==0),
                          )
    
    # add colour bar
    axes=[0.87, 0.20, 0.04, 0.6] 
    #f.subplots_adjust(top=0.95)
    f.subplots_adjust(wspace=0.001)
    f.subplots_adjust(right=axes[0]-0.03)
    cbar_ax = f.add_axes(axes)
    f.colorbar(cs, cax=cbar_ax, format=cbarform)

def vert_motion_hour(dtime=datetime(2016,1,6,7), 
        model_run='waroona_run1',
        extentname=None,
        HSkip=None):
    '''
    create vert motion slices of an hours output
    '''
    dpi=200
    if extentname is None:
        extentname=model_run.split('_')[0]
    extent=plotting._extents_[extentname]
    
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
    offset=8 if extent[0]<120 else 10
    
    for i in range(len(ff_dtimes)):
        
        utc = ff_dtimes[i]
        ltime = utc+timedelta(hours=offset)
        stitle = ltime.strftime("%Y %b %d %H:%M (LT)")
        
        # Run 4x4 and 3x3 plot creation
        for flag,subdir in zip([False,True],["4x4","3x3"]):
            vert_motion_slices(qc[i].data, w[i].data,lh,lat,lon,
                               ff=ff[i].data,
                               extentname=extentname,
                               three_by_three=flag,
                               )                
            # Save figure into animation folder with numeric identifier
            plt.suptitle(stitle)
            fio.save_fig(model_run,_sn_,utc,plt,
                         subdir=subdir,dpi=dpi)
        


### RUN THE CODE:
if __name__ == '__main__':
    # set font sizes etc.
    plotting.init_plots()
    
    mr = 'sirivan_run5_hr'
    extentname='sirivanz' # None
    hours = fio.run_info[mr]['filedates'][6:10]
    HSkip = None
    for hour in hours:
        vert_motion_hour(dtime=hour, 
                         model_run=mr, 
                         extentname=extentname, 
                         HSkip=HSkip)
    print("INFO: vert_motion_slices done")
    
    #for mr in ['waroona_old','waroona_run1']: #['sirivan_run1','waroona_old', 'waroona_run1']:
    #    hours = fio.run_info[mr]['filedates']
    #    if testing:
    #        hours = hours[:3] 
    #    for dtime in hours:
    #        vert_motion_hour(dtime, model_run=mr)
    
