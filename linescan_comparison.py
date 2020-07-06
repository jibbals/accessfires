# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:16:43 2019
    Compare model output to sirivan line scans
@author: jgreensl
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.colors as col

from glob import glob
from datetime import datetime, timedelta
from os.path import basename
from cartopy import crs as ccrs
import cartopy.io.img_tiles as cimgt

from utilities import fio, utils, plotting

_sn_ = 'linescan_comparison'

def linescan_vs_firefront(model_run='sirivan_run1'):
    """
    show linescan above firefront on as similar as possible extents
    if firefront is not avail at the desired time, show closest
    
    compares with tiff image as background, brightened as I don't know the actual wavelengths of the tiff I have...
    """
    # sirivan linescan figures extent (eyeballed using google map + matplotlib gridlines)
    linescan_extent =  plotting._extents_['sirivan_linescan']
    latlon_CRS = ccrs.PlateCarree()
        
    # read fire front (all datetimes)
    ff, sh = fio.read_fire(model_run,extent=linescan_extent,
                           firefront=True,sensibleheat=True)
    ffdates = utils.dates_from_iris(ff)
    lat = ff.coord('latitude').points
    lon = ff.coord('longitude').points
    
    # loop over linescans
    files = glob('data/linescans/SIRIVAN_*.jpg')
    files.sort()
    
    for file in files:
        fname=basename(file)
        AEST = datetime.strptime(fname,"SIRIVAN_%Y%m%d %H%M.jpg")
        # local time is offset by 10hrs from UTC
        # AEST = UTC+10H
        UTC = AEST - timedelta(hours=10)
        #dtime = local_dtime # maybe timestamp is UTC?
        dstamp = UTC.strftime('%Y%m%d %H:%M (UTC)')
        dstr = UTC.strftime('%Y%m%d_%H%M')
        
        # closest index (model output in UTC)
        di = utils.nearest_date_index(UTC,ffdates,allowed_seconds=1e10)
        ffstamp = ffdates[di].strftime('%Y%m%d %H:%M (UTC)')
        #show linescans on row 1
        img1 = mpimg.imread(file)

        fig = plt.figure(figsize=(10,13))
        ax1 = fig.add_axes([0,.5,1,.48], frameon=False)
        plt.imshow(img1)
        plt.title(dstamp)
        
        ax1.axes.get_yaxis().set_visible(False)
        ax1.axes.get_xaxis().set_visible(False)

        #show fire plan on row 2 (matching linescan times)
        # can show either googlemap (roads, rivers), or satellite map (true colour)
        # map first
        subplot_axes2 = [0,0.02,1,0.45]
        _,ax2 = plotting.map_tiff_qgis('sirivan_map_linescan.tiff', fig=fig,
                                        extent=linescan_extent,
                                        subplot_axes=subplot_axes2)
        plt.sca(ax2)
        # add firefront hourly up to current datetime indext
        # just read hourly
        hourinds = np.array([(ft.minute==0) and (ft.second==0) and ft < (ffdates[di]) for ft in ffdates])
        nt = np.sum(hourinds)
        
        # plot contours at each hour
        if nt>0:
            # contour the hours up to our linescan datetime
            for ii,dt in enumerate(ffdates[hourinds]):
                ffdata = ff[hourinds][ii].data.data
                if np.min(ffdata)>0:
                    continue
                plt.contour(lon, lat, ffdata.T, np.array([0]),
                            colors='pink', linewidths=2,alpha=0.8,
                            )

        plt.title(ffstamp)
        # final outline contour
        ffdata = ff[di].data.data
        if np.min(ffdata)<0:
            plt.contour(lon, lat, ffdata.T, np.array([0]), 
                        linestyles='dotted', alpha=0.6,
                        colors='magenta', linewidths=2)
            shdata = sh[di].data.data
            # show heat flux:
            #print("DEBUG:",np.shape(shdata),np.min(shdata),np.max(shdata))
            if np.max(shdata) > 100:
                #shdata[shdata<500] = np.NaN
                shmin,shmax=1e2,2e5
                cnorm = col.LogNorm(vmin=shmin,vmax=shmax,)
                shlevels = np.logspace(2,5.4,)
                cs=plt.contourf(lon,lat,shdata.T+1,shlevels, 
                                cmap='hot', norm=cnorm, vmin=shmin, vmax=shmax,
                                alpha=0.5)
                cax = fig.add_axes([.91,.1,.02,.3], frameon=False)
                plt.colorbar(cs,cax=cax, label='sensible heat flux',
                             ticks=[1e2,1e3,1e4,1e5,2e5],
                             orientation='vertical')
        
        fio.save_fig(model_run, _sn_, "linescan_%s"%dstr, plt)

if __name__=='__main__':
    
    linescan_vs_firefront("sirivan_run3_hr")
