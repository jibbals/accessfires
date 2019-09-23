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

def pyrocb(dtime,
           ztop=15000,
           cloud_threshold=constants.cloud_threshold,
           model_version='waroona_run1',
           ext='.png',
           dpi=300,
           ):
    '''
    311: mean vert windspeed between 500m and 1500m (?)
    312: vert motion transect
    313: vert motion transect roughly perp to 312
    INPUTS:
        datetime is used to read model and fire outputs
        extentname = { 'waroona' | 'sirivan' } choice of two fire locations
        nquivers reduces quiver density (may need to fiddle)
        quiverscale changes how long the arrows are (also may need to fiddle)
        ext is the plot extension { '.png' | '.eps' }
        dpi is plot quality (100 is default, higher is better for publication)
        old is a flag to look at pcfile output (from older model run)
    '''
    
    ### First use datetime and extentname to read correct outputs:
    extentname=model_version.split('_')[0]
    extent = plotting._extents_[extentname]

    ## read um output over extent [t, lev, lat, lon]
    cubes = fio.read_model_run(model_version, fdtime=dtime, extent=extent,
                               add_z=True, add_topog=True)
    w,  = cubes.extract('upward_air_velocity')
    qc, = cubes.extract('qc')
    topog, = cubes.extract('surface_altitude')
    z, = cubes.extract('z_th') 
    
    lat = w.coord('latitude').points
    lon = w.coord('longitude').points

    ## fire front
    ffdtimes = utils.dates_from_iris(w)
    ff1, = fio.read_fire(dtimes=ffdtimes, extent=extent, firefront=True)
    
    ff=ff1
    if model_version=="waroona_old":
        # just want 1 time step for z
        # interpolat ff onto old lats and lons    
        ff=ff1.interpolate([('longitude',w.coord('longitude').points), 
                            ('latitude',w.coord('latitude').points)],
                            iris.analysis.Linear()) 
    
    # take mean of vert motion between lvls 25-48 approx 500m - 1500m
    with warnings.catch_warnings():
        # ignore warning from collapsing non-contiguous dimension:
        warnings.simplefilter('ignore')
        wmean = w[:,25:48,:,:].collapsed('model_level_number', iris.analysis.MEAN)
    h0,h1 = wmean.coord('level_height').bounds[0]

    
    ## Plotting setup
    # set font sizes
    plotting.init_plots()
    # get plot extent, and transect
    start,end = plotting._transects_["pyrocb_waroona"]
    startx1,endx1 = plotting._transects_["pyrocbx1_waroona"]
    startx2,endx2 = plotting._transects_["pyrocbx2_waroona"]
    colorx1='b'
    colorx2='teal'
    
    # for each timestep:
    for i in range(len(ffdtimes)):
        # datetime timestamp for file,title
        stitle = ffdtimes[i].strftime("Vertical motion %Y %b %d %H:%M (UTC)")
        
        # figure setup
        f=plt.figure(figsize=[7,10])
        plt.subplot(3,1,1)
    
        ### Plot 311
        # top panel is wind speed surface values
        cs= plotting.map_contourf(extent,wmean[i].data,lat,lon,
                                  clevs = np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0])),
                                  cbar=False, 
                                  cmap=plotting._cmaps_['verticalvelocity'],
                                  norm=col.SymLogNorm(0.25),
                                  cbarform=tick.ScalarFormatter(),
                                  clabel='m/s')
        plt.title('Mean(%3.0fm - %4.0fm)'%(h0,h1))
        
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
        plotting.map_add_locations_extent(extentname)
    
        # Add fire outline
        with warnings.catch_warnings():
        # ignore warning when there are no fires:
            warnings.simplefilter('ignore')
            plt.contour(lon,lat,np.transpose(ff[i].data),np.array([0]), 
                        colors='red')
        
        ### transect plots
        ###
        ax2=plt.subplot(3,1,2)
        
        ## Plot vert motion transect
        wslice, xslice, zslice = plotting.transect_w(w[i].data,z.data,lat,lon,start,end, npoints=100,
                                                     title='',
                                                     topog=topog.data, ztop=ztop,lines=None,colorbar=False)
        ## add cloud outlines
        ## Add contour where clouds occur
        qcslice = utils.cross_section(qc[i].data,lat,lon,start,end, npoints=100)
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
        
        
        for ax3,startx,endx,colorx in zip([plt.subplot(3,2,5), plt.subplot(3,2,6)],
                                          [startx1,startx2],
                                          [endx1,endx2],
                                          [colorx1,colorx2]):
            plt.sca(ax3)
        
            ## Plot vert motion transect
            wslicex,xslicex,zslicex = plotting.transect_w(w[i].data,z.data,lat,lon,startx,endx,
                                                          title='',
                                                          npoints=100,
                                                          topog=topog.data, 
                                                          colorbar=False,
                                                          ztop=ztop,
                                                          lines=None)
            ## add cloud outlines
            ## Add contour where clouds occur
            qcslicex = utils.cross_section(qc[i].data,lat,lon,startx,endx, npoints=100)
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
        # Show transect start and end
        #xticks,xlabels = plotting.transect_ticks_labels(startx,endx)
        #plt.xticks(xticks[0::2],xlabels[0::2]) # show transect start and end
        
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
        
        # Save figure into animation folder with numeric identifier
        plt.suptitle(stitle)
        fio.save_fig(model_version, _sn_, ffdtimes[i], plt, dpi=dpi)
            
if __name__ == '__main__':
    
    print("INFO: testing cloud_outline.py")
    #emberstorm_clouds(datetime(2016,1,5,15))
    
    for dtime in [ datetime(2016,1,6,7) + timedelta(hours=x) for x in range(2) ]:
        for mv in ['waroona_run1','waroona_old']:
            pyrocb(dtime, model_version=mv)

