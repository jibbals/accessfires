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
matplotlib.use('Agg')

# plotting stuff
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.ticker as tick
import numpy as np
from datetime import datetime,timedelta
import iris # file reading and constraints etc

# local modules
from utilities import plotting, utils, fio


def emberstorm_clouds(dtime,
                      extentname='waroona',
                      vectorskip=13,
                      quiverscale=60,
                      ztop=4000,
                      ext='.png',
                      dpi=400,
                  ):
    '''
    311: mean vert windspeed between 500m and 1500m (?)
    312: vert motion transect along emberstorm
    313: vert motion transect along emberstormx (roughly perp to emberstorm)
    INPUTS:
        datetime is used to read model and fire outputs
        extentname = { 'waroona' | 'sirivan' } choice of two fire locations
        vectorskip reduces quiver density (may need to fiddle)
        quiverscale changes how long the arrows are (also may need to fiddle)
        ext is the plot extension { '.png' | '.eps' }
        dpi is plot quality (100 is default, higher is better for publication)
    '''
    
    # figure name and location
    pnames="figures/%s/emberstorm/fig_%s%s"
    
    ## time intervals for model output this hour
    ffdtimes = [ dtime + timedelta(minutes=x) for x in np.arange(10,61,10,dtype=float)]
    
    ### First use datetime and extentname to read correct outputs:
    extent = plotting._extents_[extentname]
    ## Read fire front over extent [t, lat, lon]
    ff, = fio.read_fire(ffdtimes, extent=extent, firefront=True)
    ## read um output over extent [t, lev, lat, lon]
    slv, _, th1, th2 = fio.read_waroona_iris(dtime,extent=extent)
    qc, = th2.extract('qc')
    lat = qc.coord('latitude').points
    lon = qc.coord('longitude').points
    w,  = th1.extract('upward_air_velocity')
    # take mean of vert motion between lvls 25-48 approx 500m - 1500m
    wmean = w[:,25:48,:,:].collapsed('model_level_number', iris.analysis.MEAN)
    h0,h1 = wmean.coord('level_height').bounds[0]
    #     bounds=array([[ 450.    , 1600.0004]])
    topog,= slv.extract('topog') # [ lat, lon]
    z,    = th1.extract('z_th') # [t, lev, lat, lon]
    ## Plotting setup
    # set font sizes
    plotting.init_plots()
    # get plot extent, and transect
    start,end   = plotting._transects_["emberstorm"]
    startx,endx = plotting._transects_["emberstormx"]
    
    # for each timestep:
    for i in range(6):
        # datetime timestamp for file,title
        dstamp = ffdtimes[i].strftime("%Y%m%d%H%M")
        stitle = ffdtimes[i].strftime("%Y %b %d %H:%M (UTC)")
        pname = pnames%(extentname,dstamp,ext)
        # figure setup
        plt.figure(figsize=[7,10])
        plt.subplot(3,1,1)
    
        ### Plot 311
        # top panel is wind speed surface values
        plotting.map_contourf(extent,wmean[i].data,lat,lon,
                              clevs = np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0])),
                              cmap=plotting._cmaps_['verticalvelocity'],norm=col.SymLogNorm(0.25),cbarform=tick.ScalarFormatter(),
                              clabel='m/s')
        plt.title('Vertical motion mean(%3.0fm - %4.0fm)'%(h0,h1))
        
    
        # start to end x=[lon0,lon1], y=[lat0, lat1]
        plt.plot([start[1],end[1]],[start[0],end[0], ], '--m', 
                 linewidth=2, 
                 marker='X', markersize=7)
        plt.plot([startx[1],endx[1]],[startx[0],endx[0], ], '--b', 
                 linewidth=2)
        
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

    
        # Add fire outline
        plt.contour(lon,lat,np.transpose(ff[i].data),np.array([0]), 
                    colors='red')
        
        ### transect plots
        ###
        ax2=plt.subplot(3,1,2)
        
        ## Plot vert motion transect
        plotting.transect_w(w[i].data,z.data,lat,lon,start,end,
                            topog=topog.data, ztop=ztop)
        ## add cloud outlines
        ## Add contour where clouds occur
        
        plt.title("Vertical motion (m/s)")
        plt.ylabel('height (m)')
        plt.xlabel('')
        for spine in [ax2.spines['bottom'], ax2.spines['top']]:
            spine.set_linestyle('dashed')
            spine.set_capstyle("butt")
            spine.set_color("m")
        
        ax3=plt.subplot(3,1,3)
        
        ## Plot vert motion transect
        plotting.transect_w(w[i].data,z.data,lat,lon,startx,endx,
                            topog=topog.data, ztop=ztop)
        ## add cloud outlines
        ## Add contour where clouds occur
        plt.title("Vertical motion (m/s)")
        plt.ylabel('height (m)')
        plt.xlabel('')
        for spine in [ax3.spines['bottom'], ax3.spines['top']]:
            spine.set_linestyle('dashed')
            spine.set_capstyle("butt")
            spine.set_color("blue")
        
        # Show transect start and end
        #xticks,xlabels = plotting.transect_ticks_labels(startx,endx)
        #plt.xticks(xticks[0::2],xlabels[0::2]) # show transect start and end
        
        # Save figure into animation folder with numeric identifier
        plt.suptitle(stitle)
        print("INFO: Saving figure:",pname)
        plt.savefig(pname,dpi=dpi)
        plt.close()
    
            
if __name__ == '__main__':
    
    print("INFO: testing cloud_outline.py")
    #emberstorm_clouds(datetime(2016,1,5,15))
    
    for dtime in [ datetime(2016,1,6,7) + timedelta(hours=x) for x in range(4) ]:
        emberstorm_clouds(dtime)

