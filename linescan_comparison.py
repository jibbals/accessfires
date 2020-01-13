# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:16:43 2019
    Compare model output to sirivan line scans
@author: jgreensl
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from glob import glob
from datetime import datetime, timedelta
from os.path import basename
from cartopy import crs as ccrs
import cartopy.io.img_tiles as cimgt

from utilities import fio, utils, plotting

_sn_ = 'linescan_comparison'

def linescan_vs_firefront(model_run='sirivan_run1',option='stamen'):
    """
    show linescan above firefront on as similar as possible extents
    if firefront is not avail at the desired time, show closest
    
    can show against satellite imagery or googlemap (with roads etc)
    options: 
        'google' for google map
        'satellite' for satellite map
        'topog' for model topography
        NOTE: These are currently broken... use 'topog' as option
    """
    # sirivan linescan figures extent (eyeballed using google map + matplotlib gridlines)
    linescan_extent =  [149.48, 150.04, -32.18, -31.85]
    latlon_CRS = ccrs.PlateCarree()
        
    # read fire front (all datetimes)
    ff, = fio.read_fire(model_run,extent=linescan_extent)
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
        subplot_extent2 = [0,0.02,1,0.45]
        if option=='stamen':
            stamen_terrain = cimgt.Stamen('terrain-background')
            ax2 = fig.add_axes(subplot_extent2, frameon=False, projection=stamen_terrain.crs)
            ax2.set_extent(linescan_extent)
            ax2.add_image(stamen_terrain, 10) # 10 seems to be max zoom for resolution
            mproj = stamen_terrain.crs
        else:
            _,ax2,mproj = plotting.map_tiff('sirivan', fig=fig,
                                                  extent=linescan_extent, 
                                                  subplot_row_col_n=[2,1,1])
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
                            colors='orangered', linewidths=2,
                            transform=latlon_CRS)

        # final outline contour
        ffdata = ff[di].data.data
        if np.min(ffdata)<0:
            plt.contour(lon, lat, ffdata.T, np.array([0]), 
                        linestyles='dotted',
                        colors='magenta', linewidths=2, transform=latlon_CRS)
        plt.title(ffstamp)
        fio.save_fig('sirivan_run1',_sn_,"linescan_%s"%dstr,plt,dpi=150)

if __name__=='__main__':
    
    linescan_vs_firefront()
