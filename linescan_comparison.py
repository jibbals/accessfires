# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 12:16:43 2019
    Compare model output to sirivan line scans
@author: jgreensl
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import warnings
from glob import glob
from datetime import datetime, timedelta
from os.path import basename
from cartopy import crs as ccrs

from utilities import fio, utils, plotting

_sn_ = 'linescan_comparison'

def linescan_vs_firefront(model_run='sirivan_run1',google=False):
    """
    show linescan above firefront on as similar as possible extents
    if firefront is not avail at the desired time, show closest
    
    can show against satellite imagery or googlemap (with roads etc)
    """
    # sirivan linescan figures extent (eyeballed using google map + matplotlib gridlines)
    linescan_extent =  [149.48, 150.04, -32.18, -31.85]
    uarbry_offset = timedelta(hours=10) # UARBRY LOCAL TIME OFFSET FROM UTC
    latlon_CRS = ccrs.PlateCarree()
    
    # read fire front (all datetimes)
    ff, = fio.read_fire('sirivan_run1',extent=linescan_extent)
    ffdates = utils.dates_from_iris(ff)
    lat = ff.coord('latitude').points
    lon = ff.coord('longitude').points
    
    # loop over linescans
    files = glob('data/linescans/SIRIVAN_*.jpg')
    files.sort()
    
    for file in files:
        fname=basename(file)
        local_dtime = datetime.strptime(fname,"SIRIVAN_%Y%m%d %H%M.jpg")
        # local time is offset by 10hrs from UTC
        #dtime = local_dtime - uarbry_offset
        dtime = local_dtime # maybe timestamp is UTC?
        dstamp = dtime.strftime('%Y%m%d %H:%M (UTC)')
        dstr = dtime.strftime('%Y%m%d_%H%M')
        
        # closest index:
        di = utils.nearest_date_index(dtime,ffdates,allowed_seconds=1e10)
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
        if google:
            _,ax2,mproj = plotting.map_google(linescan_extent, fig=fig, zoom=12,
                                              subplot_extent=[0,0.02,1,0.45],
                                              draw_gridlines=False)
        else:
            _,ax2,mproj,dproj = plotting.map_satellite(linescan_extent, fig=fig,
                                                       subplot_extent=[0,0.02,1,0.45])
        plt.sca(ax2)
        # add firefront hourly up to current datetime indext
        # just read hourly
        hourinds = np.array([(ft.minute==0) and (ft.second==0) and ft < (ffdates[di]) for ft in ffdates])
        nt = np.sum(hourinds)
        
        # plot contours at each hour
        with warnings.catch_warnings():
            # ignore no fire contour warning
            warnings.simplefilter('ignore')
            if nt>0:
                for ii,dt in enumerate(ffdates[hourinds]):
                    print("DEBUG:",ii,dt,ff.shape)
                    ffdata = ff[hourinds][ii].data.data
                    print("DEBUG:",np.min(ffdata))
                    plt.contour(lon, lat, ffdata.T, np.array([0]),
                                colors='orange', linewidths=1,
                                transform=latlon_CRS)
    
            # final outline contour
            plt.contour(lon,lat,ff[di].data.data.T, np.array([0]), 
                        linestyles='dotted',
                        colors='red', linewidths=2, transform=latlon_CRS)
        plt.title(ffstamp)
        fio.save_fig('sirivan_run1',_sn_,"linescan_%s"%dstr,plt,dpi=150)

if __name__=='__main__':
    print("STARTING")
    linescan_vs_firefront()