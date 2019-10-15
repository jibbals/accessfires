# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 11:09:54 2019
    Zoom in on ember storm over waroona fire
    Plot:
        311: mean vert windspeed between 500m and 1500m (?)
        312: vert motion transect along emberstorm
        313: vert motion transect along emberstormx (roughly perp to emberstorm)

@author: jgreensl
"""

import matplotlib
matplotlib.use('Agg',warn=False)

# plotting stuff
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.ticker as tick
import numpy as np
from datetime import datetime,timedelta
import iris # file reading and constraints etc
import warnings

# local modules
from utilities import plotting, utils, fio, constants

###
## GLOBALS
###
_sn_ = 'pyrocb'

# transects showing the pyroCB in model extents
__PYRO_TRANSECTS__ = {'waroona':{'X1':[[-32.89 , 115.9 ], [-32.88   , 116.19]],
                                 'X2':[[-32.95, 116.09], [-32.78   , 116.16]],
                                 'X3':[[-32.95, 116.15], [-32.78   , 116.09]]},
                      'sirivan':{'X1':plotting._transects_['sirivan1'],
                                 'X2':[[-32.1, 149.7  ], [-31.85 , 150.2 ]],
                                 'X3':[[-32.1, 150.2  ], [-31.85 , 149.7 ]]}}



def pyrocb(w, qc, z, wmean, topog, lat, lon,
           transect1,transect2, transect3,
           ff=None,
           wmeantitle='Average vertical motion',
           extentname=None,
           ztop=15000,
           cloud_threshold=constants.cloud_threshold,
           ext='.png',
           dpi=300,
           ):
    '''
    311: mean vert windspeed between 500m and 1500m (?)
    312: vert motion transect (long transect)
    325: vert motion transect (short transect 1)
    326: vert motion transect (short transect 2)
    
    INPUTS:
        w: vertical motion m/s [z, lat, lon]
        qc: cloud liquid + ice g/kg [z, lat, lon]
        z: model level heights m [z, lat, lon]
        wmean: mean vertical motion for subplot 311 m/s [lat,lon]
        topog: surface altitude m [lat,lon]
        transectN: [[lat0,lon0],[lat1,lon1]] of transect
        ff: firefront array (-ve is burnt), optional
        ztop: height of z axis in transect
        cloud_threshold: g/kg threshold of liquid + ice in air
        ext is the plot extension { '.png' | '.eps' }
        dpi is plot quality (100 is default, higher is better for publication)
    '''
    ## Plotting setup
    # set font sizes etc
    plotting.init_plots()
    # get plot extent, and transect
    extent=[lon[0],lon[-1],lat[0],lat[-1]]
    start,end = transect1
    startx1,endx1 = transect2
    startx2,endx2 = transect3
    colorx1='b'
    colorx2='teal'
    
       
    # figure setup
    f=plt.figure(figsize=[7,10])
    plt.subplot(3,1,1)

    ### Plot 311
    # top panel is wind speed surface values
    clevs_vertwind = np.union1d(np.union1d(2.0**np.arange(-2,6),
                                           -1*(2.0**np.arange(-2,6))),
                                np.array([0]))
    cs, _ = plotting.map_contourf(extent,wmean,lat,lon,
                                   clevs = clevs_vertwind,
                                   cbar=False, 
                                   cmap=plotting._cmaps_['verticalvelocity'],
                                   norm=col.SymLogNorm(0.25),
                                   cbarform=tick.ScalarFormatter(),
                                   clabel='m/s')
    plt.title(wmeantitle)
    
    # start to end x=[lon0,lon1], y=[lat0, lat1]
    plt.plot([start[1],end[1]],[start[0],end[0], ], '--m', 
             linewidth=2)
    def myarrow(yx0, yx1):
        y0,x0=yx0
        y1,x1=yx1
        dx,dy = (x1-x0)/15.0,(y1-y0)/15.0
        plt.arrow(x0-1.4*dx,y0-1.4*dy,dx,dy,width=0.0015)
    myarrow(start,end)
    for startx,endx,xcolor in zip([startx1,startx2],[endx1,endx2],[colorx1,colorx2]):
        plt.plot([startx[1],endx[1]],[startx[0],endx[0], ], '--',
                 color=xcolor,
                 linewidth=2)
        myarrow(startx,endx)
    
    # add nearby towns
    if extentname is not None:
        plotting.map_add_locations_extent(extentname)

    # Add fire outline
    if ff is not None:
        with warnings.catch_warnings():
        # ignore warning when there are no fires:
            warnings.simplefilter('ignore')
            plt.contour(lon,lat,np.transpose(ff),np.array([0]), colors='red')
    
    ### transect plots
    ###
    ax2=plt.subplot(3,1,2)
    
    ## Plot vert motion transect
    wslice, xslice, zslice = plotting.transect_w(w, z, lat, lon, start, end,
                                                 npoints=100, title='',
                                                 topog=topog, ztop=ztop,
                                                 lines=None, colorbar=False)
    ## add cloud outlines
    ## Add contour where clouds occur
    qcslice = utils.cross_section(qc, lat, lon, start, end, npoints=100)
    with warnings.catch_warnings():
        # ignore warning when there are no clouds:
        warnings.simplefilter('ignore')
        plt.contour(xslice,zslice,qcslice,np.array([cloud_threshold]),colors='k')

    plt.ylabel('height (m)')
    plt.xlabel('')
    for spine in [ax2.spines['bottom'], ax2.spines['top']]:
        spine.set_linestyle('dashed')
        spine.set_capstyle("butt")
        spine.set_color("m")
        spine.set_linewidth(2)
    
    ax31, ax32 = [plt.subplot(3,2,5), plt.subplot(3,2,6)]
    for ax3,startx,endx,colorx in zip([ax31, ax32],
                                      [startx1,startx2],
                                      [endx1,endx2],
                                      [colorx1,colorx2]):
        plt.sca(ax3)
    
        ## Plot vert motion transect
        wslicex,xslicex,zslicex = plotting.transect_w(w, z, lat, lon, 
                                                      startx, endx,
                                                      title='',
                                                      npoints=100,
                                                      topog=topog, 
                                                      colorbar=False,
                                                      ztop=ztop,
                                                      lines=None)
        ## add cloud outlines
        ## Add contour where clouds occur
        qcslicex = utils.cross_section(qc, lat, lon, startx, endx, 
                                       npoints=100)
        with warnings.catch_warnings():
            # ignore warning when there are no clouds:
            warnings.simplefilter('ignore')
            plt.contour(xslicex,zslicex,qcslicex,np.array([cloud_threshold]),colors='k')

        plt.xlabel('')
        for spine in [ax3.spines['bottom'], ax3.spines['top']]:
            spine.set_linestyle('dashed')
            spine.set_capstyle("butt")
            spine.set_color(colorx)
            spine.set_linewidth(2)
    
    # remove y ticks from last figure
    plt.yticks([],[])
        
    # add colour bar
    axes=[0.85, 0.20, 0.03, 0.6] 
    #f.subplots_adjust(top=0.95)
    #f.subplots_adjust(hspace=0.05)
    f.subplots_adjust(right=axes[0]-0.01)
    cbar_ax = f.add_axes(axes)
    cb = f.colorbar(cs, cax=cbar_ax, 
                    format=tick.ScalarFormatter())
    # -ve labelpad moves label closer to colourbar
    cb.set_label('m/s', labelpad=-3)
        
            
def pyrocb_model_run(model_run='waroona_run1', dtime=datetime(2016,1,5,15)):
    """
    """
    ### First use datetime and extentname to read correct outputs:
    extentname=model_run.split('_')[0]
    extent = plotting._extents_[extentname]
    X1 = __PYRO_TRANSECTS__[extentname]['X1']
    X2 = __PYRO_TRANSECTS__[extentname]['X2']
    X3 = __PYRO_TRANSECTS__[extentname]['X3']
    
    ## read um output over extent [t, lev, lat, lon]
    cubes = fio.read_model_run(model_run, fdtime=[dtime], extent=extent,
                               add_z=True, add_topog=True)
    
    w,  = cubes.extract('upward_air_velocity')
    qc, = cubes.extract('qc')
    topog, = cubes.extract('surface_altitude')
    z, = cubes.extract('z_th') 
    
    lat = w.coord('latitude').points
    lon = w.coord('longitude').points

    ## fire front
    ffdtimes = utils.dates_from_iris(w)
    ff, = fio.read_fire(model_run=model_run, dtimes=ffdtimes, extent=extent, firefront=True)
    
    
    # take mean of vert motion between lvls 25-48 approx 500m - 1500m
    with warnings.catch_warnings():
        # ignore warning from collapsing non-contiguous dimension:
        warnings.simplefilter('ignore')
        wmean = w[:,25:48,:,:].collapsed('model_level_number', iris.analysis.MEAN)
    h0,h1 = wmean.coord('level_height').bounds[0]
    
    # for each timestep:
    for i in range(len(ffdtimes)):
        # datetime timestamp for file,title
        stitle = ffdtimes[i].strftime("Vertical motion %Y %b %d %H:%M (UTC)")
        wmeantitle='Mean(%3.0fm - %4.0fm)'%(h0,h1)
        
        #print("DEBUG: =====")
        #print(cubes)
        #print(ff.shape)
        #print(wmean.shape)
        #print(i)
        #print("DEBUG: =====")
        pyrocb(w=w[i].data.data, qc=qc[i].data.data, z=z.data.data, 
               wmean=wmean[i].data.data, topog=topog.data.data,
               lat=lat,lon=lon, transect1=X1, transect2=X2, transect3=X3,
               ff=ff[i].data.data,
               wmeantitle=wmeantitle,
               extentname=extentname)
        # Save figure into animation folder with numeric identifier
        plt.suptitle(stitle)
        fio.save_fig(model_run, _sn_, ffdtimes[i], plt, dpi=200)

if __name__ == '__main__':
    
    testing=False
    model_runs = ['sirivan_run1', 'waroona_run1','waroona_old']
    if testing:
        model_runs=model_runs[0:2]
        
    for mr in model_runs :
        dtimes = fio.model_outputs[mr]['filedates']
        if testing:
            dtimes = dtimes[0:2]
        for dtime in dtimes:
            pyrocb_model_run(model_run=mr, dtime=dtime)
