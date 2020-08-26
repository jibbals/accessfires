#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 15:45:40 2019
    
    Compare two model runs together...
    
@author: jesse
"""

import matplotlib
matplotlib.use("Agg",warn=False)

from matplotlib import colors, ticker, patheffects

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from datetime import datetime,timedelta
from scipy.stats import gaussian_kde, cumfreq
import warnings # silence some contour warnings

from utilities import utils, plotting, fio, constants
import fireplan
import emberstorm

## GLOBAL
#Script name
_sn_ = 'run_comparison'

def compare_surface(mr1='waroona_run2', mr2='waroona_run2uc', hour=datetime(2016,1,5,15)):
    """
        3 by 3 plot, ROWS: H winds, V winds, Theta, COLS: Run2, Run2UC, DIFF
    """
    fig, axes = plt.subplots(3,3)

def compare_winds(mr1='waroona_run2', mr2='waroona_run2uc', 
                  hour=datetime(2016,1,5,15),
                  extent=None,
                  subsubdir=None):
    """
        Compare winds between two different runs, looking at vertically binned averages
        first row shows model run 1 (and figures are saved under this run)
        second row shows model run 2
        final row(s) show density or cumulative density from averaged metric
        
        two separate folders:
            horizontal/
                4 averaged binds (altitudes 0-500, 500-2000, 2000-5000, 5000-9000)
                mr1
                mr2
                distribution histograms for wind speed and direction
            vertical
                same as horizontal, with cdf for wind speed
    """
    extentname = mr1.split('_')[0]
    if extent is None:
        extent = plotting._extents_[extentname]
    

    cubes1 = fio.read_model_run(mr1, hour, extent=extent, add_winds=True)
    cubes2 = fio.read_model_run(mr2, hour, extent=extent, add_winds=True)
    
    # pull out horizontal winds
    u,v,s,wd = cubes1.extract(['u','v','s','wind_direction'])
    cu,cv,cs,cwd = cubes2.extract(['u','v','s','wind_direction'])
    dates = utils.dates_from_iris(u)
    height = utils.height_from_iris(s)#.coord('level_height').points
    lats = s.coord('latitude').points
    lons = s.coord('longitude').points
    clats = cs.coords('latitude').points
    clons = cs.coords('longitude').points

    ff1,ff2=None,None
    if fio.run_info[mr1]['hasfire']:
        ff1, = fio.read_fire(mr1,dtimes=dates,extent=extent,firefront=True)
    if fio.run_info[mr2]['hasfire']:
        ff2, = fio.read_fire(mr2,dtimes=dates,extent=extent,firefront=True)
    
    # density plot arguments
    bandwidth=1
    # x axis for wind speed
    xs = np.linspace(0,20,100)
    # x axis for wind dir
    xwd = np.linspace(0,360,100)
    # x axis for vert wind speed
    xw = np.linspace(-32,32,100)
    
    ## Colourmap setup for contourf plots
    # Horizontal wind colorbar should be the same for both runs
    hmaxthresh=[5,10,15,20,25,30,35]
    hcmap=plotting._cmaps_['windspeed']
    # vertical wind colourbar is constant
    wcmap=plotting._cmaps_['verticalvelocity']
    wnorm=colors.SymLogNorm(0.25) # linear to +- 0.25, then log scale
    wcontours=np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0]))
    
    # make 4 vertical bins
    row1 = (0<=height) * (height<500)
    row2 = (500<=height) * (height<2000)
    row3 = (2000<=height) * (height<5000)
    row4 = (5000<=height) * (height<9000)
    
    # loop over datetimes:
    for di, date in enumerate(dates):
        
        s1, s2, s3, s4 = [np.mean(s[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        u1, u2, u3, u4 = [np.mean(u[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        v1, v2, v3, v4 = [np.mean(v[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        wd1, wd2, wd3, wd4 = [np.mean(wd[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        
        cs1, cs2, cs3, cs4 = [np.mean(cs[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        cu1, cu2, cu3, cu4 = [np.mean(cu[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        cv1, cv2, cv3, cv4 = [np.mean(cv[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        cwd1, cwd2, cwd3, cwd4 = [np.mean(cwd[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        
        
        # Plotting
        plt.close()
        fig, axes = plt.subplots(4,4,figsize=[12,12])
        for i, (si, ui, vi, wdi, csi, cui, cvi, cwdi) in \
            enumerate(zip([s1,s2,s3,s4],[u1,u2,u3,u4],[v1,v2,v3,v4],[wd1,wd2,wd3,wd4],
                  [cs1,cs2,cs3,cs4],[cu1,cu2,cu3,cu4],[cv1,cv2,cv3,cv4],[cwd1,cwd2,cwd3,cwd4])):
            ## Show contourf of wind speed
            plt.sca(axes[0,i])
            
            # first determine colourmap contours so that they are the same between plots (and useful)
            hmax_index=np.sum(np.max(np.max(si),np.max(csi))>hmaxthresh)
            hcontours=np.linspace(0,hmaxthresh[hmax_index],20)
            
            # plot the filled contour for h-wind speeds
            plotting.map_contourf(extent, si, lats, lons, cmap=hcmap,
                                  clabel="", cbar=False,
                                  cbarform=None,
                                  levels=hcontours)
            # overlaid with quiver of wind dir
            #plotting.map_quiver(ui,vi,lats,lons,nquivers=7)
            plt.streamplot(lons,lats,ui,vi,color='grey',minlength=0.5)
            # set limits back to latlon limits
            plt.ylim(lats[0],lats[-1])
            plt.xlim(lons[0],lons[-1])
            
            # Add locations and fire
            plotting.map_add_locations_extent(extentname, 
                                              hide_text=True,
                                              color='k')
            if ff1 is not None:
                plotting.map_fire(ff1[di].data,lats,lons)
            # Add ylabel on left most plot
            if i==0: plt.ylabel(mr1)
            # add title along top row
            plt.title(['<500m','500m-2000m','2km-5km','5km-9km'][i])
            
            ## for comparison model also
            plt.sca(axes[1,i])
            img,_ = plotting.map_contourf(extent, csi, clats, clons, cmap=hcmap, 
                                          clabel="",  
                                          cbar=False, cbarform=None,
                                          levels=hcontours,)
            
            # overlaid with quiver of wind dir
            #plotting.map_quiver(cui,cvi,lats,lons,nquivers=7)
            
            plt.streamplot(clons,clats,cui,cvi,color='grey',minlength=0.5)
            # set limits back to latlon limits
            plt.ylim(lats[0],lats[-1])
            plt.xlim(lons[0],lons[-1])
            
            plotting.map_add_locations_extent(extentname, 
                                              hide_text=True,
                                              color='k')
            # add fire
            if ff2 is not None:
                plotting.map_fire(ff2[di].data,clats,clons)
            # add label on leftmost
            if i==0: plt.ylabel(mr2)
            
            ## Add colourbar for column
            cbar_ax = fig.add_axes([0.14+0.2025*i, 0.508, 0.14, 0.01])# X Y Width Height
            cbar=fig.colorbar(img, cax=cbar_ax, format=ticker.ScalarFormatter(), pad=0,
                              orientation='horizontal')
            cbar.set_ticks(np.arange(0,hmaxthresh[hmax_index]+1,5))
            cbar.set_ticklabels(np.arange(0,hmaxthresh[hmax_index]+1,5))
            
            ## add density plot for wind speed
            plt.sca(axes[2,i])
            # density function:
            sdens = gaussian_kde(si.flatten(),bw_method=bandwidth)
            plt.plot(xs,sdens(xs), label=mr1,linewidth=2)
            csdens = gaussian_kde(csi.flatten(),bw_method=bandwidth)
            plt.plot(xs,csdens(xs), label=mr2)
            plt.yticks([],[])
            if i==0: plt.ylabel('wind speed density')
            if i==3: plt.legend()
            
            
            ## add density plot for wind dir
            plt.sca(axes[3,i])
            # density function:
            wddens = gaussian_kde(wdi.flatten(),bw_method=bandwidth)
            plt.plot(xwd,wddens(xwd), label=mr1, linewidth=2)
            cwddens = gaussian_kde(cwdi.flatten(),bw_method=bandwidth)
            plt.plot(xwd,cwddens(xwd), label=mr2, linewidth=2, linestyle='--')
            plt.yticks([],[])
            if i==0: plt.ylabel('wind dir density')
            if i==3: plt.legend()
            
        plt.suptitle(date.strftime("%Y%m%d %H:%M(UTC)"))
        subdir='horizontal'
        if subsubdir is not None:
            subdir += '/'+subsubdir
        fio.save_fig(mr2, _sn_, date, plt, subdir=subdir)
    
        ## Also want to look at vertical winds
        w, = cubes1.extract(['upward_air_velocity'])
        cw, = cubes2.extract(['upward_air_velocity'])
        
        w1, w2, w3, w4 = [np.mean(w[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        cw1, cw2, cw3, cw4 = [np.mean(cw[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        
        # Plotting
        plt.close()
        fig, axes = plt.subplots(3,4,figsize=[12,12])
        for i, (wi, cwi) in enumerate(zip([w1,w2,w3,w4],[cw1,cw2,cw3,cw4])):
            ## Show contourf of wind speed
            plt.sca(axes[0,i])
            plotting.map_contourf(extent, wi, lats,lons,cmap=wcmap,clabel="",levels=wcontours,norm=wnorm,cbar=False,cbarform=None)
            plotting.map_add_locations_extent(extentname, 
                                              hide_text=True,
                                              color='k')
            if ff1 is not None:
                plotting.map_fire(ff1[di].data,lats,lons)
            if i==0: plt.ylabel(mr1)
            plt.title(['<500m','500m-2000m','2km-5km','5km-9km'][i])
            
            ## for comparison model also
            plt.sca(axes[1,i])
            img,_=plotting.map_contourf(extent, cwi, clats,clons,cmap=wcmap,clabel="",levels=wcontours,norm=wnorm,cbar=False,cbarform=None)
            plotting.map_add_locations_extent(extentname, 
                                              hide_text=True,
                                              color='k')
            if ff2 is not None:
                plotting.map_fire(ff2[di].data,clats,clons)
            if i==0: plt.ylabel(mr2)
            
            ## add density plot for wind speed
            plt.sca(axes[2,i])
            ## density function:
            #wdens = gaussian_kde(wi.flatten(),bw_method=bandwidth)
            #plt.plot(xw,wdens(xw), label=mr1, linewidth=2)
            #cwdens = gaussian_kde(cwi.flatten(),bw_method=bandwidth)
            #plt.plot(xw,cwdens(xw), label=mr2, linewidth=2, linestyle='--')
            ## cumulative frequency
            wcdf = cumfreq(wi.flatten(), numbins=len(xw), defaultreallimits=(min(xw),max(xw)))
            plt.plot(xw,wcdf.cumcount,label=mr1, linewidth=2)
            cwcdf = cumfreq(cwi.flatten(), numbins=len(xw), defaultreallimits=(min(xw),max(xw)))
            plt.plot(xw,cwcdf.cumcount,label=mr2, linewidth=2, linestyle='--')
            plt.yticks([],[])
            if i==0: plt.ylabel('wind speed cdf')
            if i==3: plt.legend()
        ## Add colourbar
        cbar_ax = fig.add_axes([.35, 0.37, 0.3, 0.015])# X Y Width Height
        cbar=fig.colorbar(img, cax=cbar_ax, format=ticker.ScalarFormatter(), pad=0,
                          orientation='horizontal')
        cbar.set_ticks([-32,-8,-2,0,2,8,32])
        cbar.set_ticklabels([-32,-8,-2,0,2,8,32])
        # title and save
        plt.suptitle(date.strftime("%Y%m%d %H:%M(UTC)"))
        
        subdir='vertical'
        if subsubdir is not None:
            subdir += '/'+subsubdir
        fio.save_fig(mr2, _sn_, date, plt, subdir=subdir)
    
def compare_clouds(mr1='waroona_run2', mr2='waroona_run2uc', 
                   hour=datetime(2016,1,5,15), 
                   cloud_threshold=constants.cloud_threshold,
                   extent=None,
                   subsubdir=None):
    """
    Plot top-down view of cloud content within summed vertical bins
        also show distributions
    """
    # x axis for qc-max density
    colormax=0.5
    xs = np.linspace(0,colormax,200)
    ## Colourmap setup for contourf plots
    # linear between 0 and 0.01
    # log between 0.01 and 0.3
    logmax=np.log(colormax)/np.log(10)-0.01
    cmap=plotting._cmaps_['qc']
    norm=colors.SymLogNorm(0.01, vmin=0,vmax=colormax)
    clevs=np.union1d(np.union1d(np.logspace(-2,logmax,30),0),colormax) 
    
    extentname = mr1.split('_')[0]
    if extent is None:
        extent = plotting._extents_[extentname]
    
    ## Read a model run
    cubes1 = fio.read_model_run(mr1, hour, extent=extent)
    cubes2 = fio.read_model_run(mr2, hour, extent=extent)
    
    # pull out clouds
    qc, = cubes1.extract(['qc'])
    cqc, = cubes2.extract(['qc'])
    dates = utils.dates_from_iris(qc)
    height = utils.height_from_iris(qc)
    lats = qc.coord('latitude').points
    lons = qc.coord('longitude').points
    clats = cqc.coord('latitude').points
    clons = cqc.coord('longitude').points
    ff1, = fio.read_fire(mr1,dtimes=dates,extent=extent,firefront=True)
    ff2, = fio.read_fire(mr2,dtimes=dates,extent=extent,firefront=True)
    
    # density plot bandwidth
    bandwidth=1
    
    # make 4 vertical bins
    row1 = (2000<=height) * (height<3000)
    row2 = (3000<=height) * (height<5000)
    row3 = (5000<=height) * (height<8000)
    row4 = (8000<=height) * (height<15000)
    titles = ['2km-3km','3km-5km','5km-8km','8km-15km']
    
    # loop over datetimes:
    for di, date in enumerate(dates):
        
        qc1, qc2, qc3, qc4 = [np.sum(qc[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        cqc1, cqc2, cqc3, cqc4 = [np.sum(cqc[di,row,:,:].data, axis=0) for row in [row1,row2,row3,row4]]
        
        # Plotting
        plt.close()
        fig, axes = plt.subplots(3,4,figsize=[12,11])
        for i, (qci, cqci) in enumerate(zip([qc1, qc2, qc3, qc4],[cqc1, cqc2, cqc3, cqc4])):
            ## Show contourf of cloud
            plt.sca(axes[0,i])
            
            plotting.map_contourf(extent, qci, lats, lons, cmap=cmap, 
                                  clabel="", norm=norm, cbar=False, levels=clevs,
                                  cbarform=None, extend='max')
            
            # overlaid with cloud thresh line
            if np.max(qci)>cloud_threshold:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    plt.contour(lons,lats, qci, np.array([cloud_threshold]),
                                colors='teal', linewidths=2)
            
            # Add locations and fire
            plotting.map_add_locations_extent(extentname, hide_text=True)
            if ff1 is not None:
                plotting.map_fire(ff1[di].data,lats,lons)
            # Add ylabel on left most plot
            if i==0: plt.ylabel(mr1)
            # add title along top row
            plt.title(titles[i])
            
            ## for comparison model also
            plt.sca(axes[1,i])
            img,_ = plotting.map_contourf(extent, cqci, clats, clons, cmap=cmap, 
                                          norm=norm, clabel="", levels=clevs, 
                                          cbar=False, cbarform=None, extend='max')
            
            # overlaid with cloud thresh line
            if np.max(cqci)>cloud_threshold:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    plt.contour(clons, clats, cqci, np.array([cloud_threshold]),
                                colors='teal', linewidths=2)
            # add fire
            if ff2 is not None:
                plotting.map_fire(ff2[di].data,clats,clons)
            # add label on leftmost
            if i==0: plt.ylabel(mr2)
            
            ## add density plot
            plt.sca(axes[2,i])
            
            # density function:
            # only if non zero cloud content
            flag=2
            if not np.isclose(np.max(qci),0):
                qcdens = gaussian_kde(qci.flatten(), bw_method=bandwidth)
                plt.plot(xs,qcdens(xs),linewidth=2,color='k')
                flag-=1
            if not np.isclose(np.max(cqci),0):
                cqcdens = gaussian_kde(cqci.flatten(), bw_method=bandwidth)
                plt.plot(xs,cqcdens(xs),color='r', linestyle='--')
                flag-=1
            plt.yticks([],[])
            if i==0: plt.ylabel('cloud density')
            if i==3: 
                # manually add legend
                lines = [Line2D([0], [0], color='k', linewidth=2, linestyle='-'),
                         Line2D([0], [0], color='r', linewidth=1, linestyle='--')]
                labels = [mr1,mr2]
                plt.legend(lines, labels)

            
        ## Add colourbar
        cbar_ax = fig.add_axes([0.35, 0.367, 0.31, 0.01])# X Y Width Height
        cbar=fig.colorbar(img, cax=cbar_ax, format=ticker.ScalarFormatter(), 
                          pad=0, orientation='horizontal')
        cbar.set_ticks([0,0.01,0.05,0.1,0.2,0.4])
        cbar.set_ticklabels([0,0.01,0.05,0.1,0.2,0.4])
        plt.suptitle(date.strftime("Max cloud content:  %Y%m%d %H:%M(UTC)"),
                     fontsize=20)
        subdir='clouds'
        if subsubdir is not None:
            subdir += '/'+subsubdir
        fio.save_fig(mr2, _sn_, date, plt, subdir=subdir)

def compare_at_site(mr1='waroona_run2', mr2='waroona_run2uc', latlon = plotting._latlons_['AWS_wagerup']):
    """
    look at time series of metrics from two runs, compared to AWS if possible
    """
    print("TBD")
        
def compare_fire_spread(mrlist, mrcolors=None, extent=None,
                        overplot=False, HSkip=None):
    '''
        Show hourly fire spread overplotted or sidebyside
    '''
    extentname = mrlist[0].split('_')[0]
    if extent is None:
        extent = plotting._extents_[extentname]
    
    ## Plot fireplan for high res run
    # read all the fire data
    if not overplot:
        fig = plt.figure(figsize=(10,14))
        for i,mr in enumerate(mrlist):
            ff, = fio.read_fire(model_run=mr, dtimes=None,
                                extent=extent, firefront=True,
                                HSkip=HSkip)
            # ff is [time,lon,lat]
            # find most eastern part for each hour
            skip=30 # actually do 30 mins
            ffmin = np.min(ff.data[::skip],axis=1) # remove lat dimension
            fflons = ff.coord('longitude').points
            fftimes = utils.dates_from_iris(ff)[::skip]
            print("INFO: ", mr, ff.shape)
            print("INFO:    MIN LON      |     MAX LON ")
            for ti in range(ffmin.shape[0]):# every x mins
                if np.sum(ffmin[ti] < 0) > 0:
                    minlon = fflons[ffmin[ti] < 0][0]
                    maxlon = fflons[ffmin[ti] < 0][-1]
                    print("%02d%02d     %.3f      |     %.3f   "%(fftimes[ti].hour,fftimes[ti].minute,minlon,maxlon)) 


            subplot_row_col_n=[len(mrlist),1,i+1]
            _,ax = fireplan.fireplan(ff, show_cbar=False, 
                                     fig=fig, 
                                     subplot_row_col_n=subplot_row_col_n,)
            plt.title(mr)
        fio.save_fig('comparison','compare_fire_spread', 'firespread.png', plt)


def compare_transects(run1, run2, hours=[12,], extent=None, ztop=1000,
                      columntitles=['coupled','uncoupled'],
                      subsubdir=None):
    """
    examine some transects side by side for two runs
    PLOTS: 2 columns: one for each run, 4 rows:
        top row: topdown topog+winds+fire+3 transect lines
        3 more rows: transect views of wind (emberstorm_transect method)
    Transect lines roughly spread based on extent middle + 50% of the distance to extent edges
    """
    if extent is None:
        extent = plotting._extents_[run1.split('_')[0]]
    
    # transects = extent centre +- frac of extent width/height
    y0,x0 = (extent[2]+extent[3])/2.0, (extent[0]+extent[1])/2.0
    # do middle point +- 1/3 of distance to extent edge
    dx = (extent[1]-extent[0]) / 3.0
    dy = (extent[3]-extent[2]) / 3.0
    # transects: [[[lat,lon],[lat,lon]],...]
    transects = [[[y0+dy,x0-dx], [y0+dy, x0+dx]],
                 [[y0,x0-dx], [y0, x0+dx]],
                 [[y0-dy,x0-dx], [y0-dy, x0+dx]],
                 ]
    tcolors = ['teal','blue','magenta']
    
    dtimes=fio.run_info[run1]['filedates'][np.array(hours)]
    
    ## for each hour
    for dti, dt in enumerate(dtimes):
        ## read hour of data for both runs
        cubes1 = fio.read_model_run(run1, fdtime=[dt], extent=extent, 
                                   add_topog=True, add_winds=True,
                                   add_z=True, add_theta=True)
                
        u1,v1,w1,z1 = cubes1.extract(["u","v","upward_air_velocity","z_th"])
        ctimes=utils.dates_from_iris(w1)
        zd1 = z1.data.data
        theta1, = cubes1.extract("potential_temperature")
        topog1=cubes1.extract("surface_altitude")[0].data
        
        # read fire outputs
        ff1,sh1,u101,v101 = fio.read_fire(model_run=run1,
                                      dtimes=ctimes, 
                                      extent=extent,
                                      sensibleheat=True,
                                      wind=True)
        
        cubes2 = fio.read_model_run(run2, fdtime=[dt], extent=extent, 
                                   add_topog=True, add_winds=True,
                                   add_z=True, add_theta=True)
        u2,v2,w2,z2 = cubes2.extract(["u","v","upward_air_velocity","z_th"])
        zd2 = z2.data.data
        theta2, = cubes2.extract("potential_temperature")
        topog2=cubes2.extract("surface_altitude")[0].data
        # read fire outputs
        ff2,sh2,u102,v102 = fio.read_fire(model_run=run2,
                                      dtimes=ctimes, 
                                      extent=extent,
                                      sensibleheat=True,
                                      wind=True)
        
        
        
        for cti,ct in enumerate(ctimes):
        
            ## FIGURE BEGIN:
            fig = plt.figure(figsize=[13,14])
            ## looking at two separate runs
            for runi, [sh,ff,u10,v10,u,v,w,zd,theta,topog] in enumerate(zip([sh1,sh2],[ff1,ff2],[u101,u102],[v101,v102],[u1,u2],[v1,v2],[w1,w2],[zd1,zd2],[theta1,theta2],[topog1,topog2])):
                # read datacubes

                lats = ff.coord('latitude').points
                lons = ff.coord('longitude').points
                
                print("DEBUG:",cti, ct)
                print("     :",sh)
                print("     :",w)
                shd = sh[cti].data.data
                LT = ct + timedelta(hours=8)
                
                ffd = ff[cti].data.data
                u10d = u10[cti].data.data
                v10d = v10[cti].data.data
                
                
                #plt.subplot(4,2,1+runi) # top row
                _,ax = emberstorm.topdown_emberstorm(
                    fig=fig, subplot_row_col_n=(4,2,1+runi),
                    extent=extent, lats=lats, lons=lons,
                    ff=ffd, sh=shd, u10=u10d, v10=v10d,
                    topog=topog,
                    annotate=False,
                    sh_kwargs={'colorbar':False}
                    )
                if runi==1:
                    plt.yticks([],[])
                plt.title(columntitles[runi])
                for transect, tcolor in zip(transects,tcolors):
                    ## Add dashed line to show where transect will be
                    start,end = transect
                    # outline the dashed line
                    line_effects=[patheffects.Stroke(linewidth=4, foreground='darkgrey'), patheffects.Normal()]
                    plt.plot([start[1],end[1]],[start[0],end[0], ], 
                            linestyle='--',
                            color=tcolor,
                            linewidth=2, 
                            path_effects=line_effects,
                            alpha=0.9)
                
                ## Plot title
                plt.suptitle(LT.strftime('%b %d, %H%M (UTC+8)'))
                
                # Transect plots
                for trani,transect in enumerate(transects):
                    ax=plt.subplot(4,2,3+runi+2*trani)
                    
                    trets = emberstorm.transect_emberstorm(
                        u[cti].data.data,
                        v[cti].data.data,
                        w[cti].data.data,
                        zd, lats, lons, transect,
                        topog=topog,
                        sh=shd,
                        theta=theta[cti].data.data,
                        ztop=ztop,
                        theta_contourargs={'levels':[297,300,303,306,309]},
                        wind_contourargs={'colorbar':False},
                        )
                    if runi==1:
                        plt.yticks([],[])
                    
                    # match splines to transect lines
                    plotting.set_spine_color(ax,tcolors[trani])
                    # finally add desired annotations
                    plotting.annotate_max_winds(trets['s'])
                    
            
            ## add wind speed colorbar
            cbar_ax = fig.add_axes([0.92, 0.4, 0.02, 0.2]) # X Y Width Height
            cmap = matplotlib.cm.get_cmap(name=plotting._cmaps_['windspeed'])
            #cmap.set_over('1.0')
            #norm = matplotlib.colors.Normalize(vmin=0, vmax=22.5)
            norm = matplotlib.colors.BoundaryNorm(np.arange(0,22.51,2.5), cmap.N)
            
            cb1 = matplotlib.colorbar.ColorbarBase(cbar_ax, cmap=cmap,
                                                   norm=norm,
                                                   extend='max',
                                                   orientation='vertical')
            cb1.set_label('ms$^{-1}$')
            
                
            ## SAVE FIGURE
            subdir='transects'
            if subsubdir is not None:
                subdir += '/'+subsubdir
            plt.tight_layout(rect=[0,0.03,0.92,0.975]) #(left, bottom, right, top)
            fio.save_fig(run2,_sn_,ct,subdir=subdir,plt=plt)
    

if __name__=='__main__':
    
    ext_pcb=plotting._extents_['waroonaz']
    ext_es1=emberstorm._emberstorm_centres_['waroona_run3']['first']['extent']
    ext_es2=emberstorm._emberstorm_centres_['waroona_run3']['second']['extent']
    hrs_pcb=range(9,19)
    hrs_es1=range(16,24)
    hrs_es2=range(26,35)
    ## Compare transects
    if False:
        # Look at some different spots for es2
        ext_es2s=[[115.7, 116.0, -33.02, -32.84], # nearly like std one
                [115.72,115.95, -33.0, -32.87], # zoomed in
                [115.72,115.91, -33.0, -32.9], # zoomed in further
                ]
        ssd_es2s=['test1','test2','test3']
        #for extent, hours, subsubdir in zip([ext_pcb,ext_es1,ext_es2],[hrs_pcb,hrs_es1,hrs_es2],['pcb','es1','es2']):
        for extent, hours, subsubdir in zip(ext_es2s,[hrs_es2,hrs_es2,hrs_es2],ssd_es2s):
            compare_transects('waroona_run3','waroona_run3uc', 
                              extent=extent,
                              hours=hours, 
                              subsubdir=subsubdir,
                              columntitles=['coupled','uncoupled'],
                              )
        # look at emberstorm area
    
    ## Compare firefronts
    ## TODO NEEDS UPDATING to fix contour colouring
    if False:
        #compare_fire_spread(['sirivan_run1','sirivan_run1_hr', 'sirivan_run3_hr'], HSkip=7)
        compare_fire_spread(['waroona_run3','waroona_run3e'], HSkip=7)
    
    ## Compare topdown views of winds and clouds
    if False:
        mr1,mr2 = ['waroona_run3','waroona_run3uc']
        extent=emberstorm._emberstorm_centres_[mr1]['second']['extent']
        #hours = emberstorm._emberstorm_centres_[mr1]['first']['hours']
        hours = range(40,48)
        subsubdir='emberstorm2'
        for fdate in fio.run_info[mr1]['filedates'][hours]:
            compare_winds(mr1=mr1,mr2=mr2, hour=fdate,
                          extent=extent,
                          subsubdir=subsubdir)
            compare_clouds(mr1=mr1,mr2=mr2, hour=fdate,
                           extent=extent,
                           subsubdir=subsubdir)
        
    ## look at overall burn area
    # Lets loop over and compare run3 with it's uncoupled brother
    if True:
        mr1,mr2 = ['sirivan_run4','sirivan_run5_hr']
        extent=None # default extent applied
        extent=[149.6,150,-32.1,-31.9]
        hrs = range(8,23)
        compare_transects(mr1,mr2, hours=hrs, columntitles=[mr1,mr2],extent=extent)
        for di,fdate in enumerate(fio.run_info[mr1]['filedates']):
            if fdate not in fio.run_info[mr2]['filedates']:
                print("INFO: skipping %s: it is in %s but not in %s"%(fdate.strftime("%Y%m%d"),mr1,mr2))
                continue
            compare_winds(mr1=mr1,mr2=mr2, hour=fdate,extent=extent)
            compare_clouds(mr1=mr1,mr2=mr2, hour=fdate,extent=extent)
    
    print("INFO: run_comparison.py done")
