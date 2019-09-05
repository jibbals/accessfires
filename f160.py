# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 11:18:23 2019
    Create F160 or skewt logp plots for model output
    read pressure and temperature, 
    estimate boundary level stability,height,etc
    estimate fire parameters for something(?)
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


def f160(press,temps, tempd, latlon, nearby=2, alpha=0.4):
    '''
    show skewt logp plot of temperature profile at latlon
    Does model give us dewpoint or equiv to also plot here??
    input cubes: pressures, temperatures are [z,lat,lon]
    profile for cubes interpolated to latlon point will be shown with linewidth of 2
    nearest values to lat lon within <nearby> indices will be plotted at low alpha
    '''
    
    ## first interpolate pressure and temperature to latlon
    lons=press.coord('longitude')
    lats=press.coord('latitude')
    press.convert_units('hPa') # convert to hPa
    temps.convert_units('K') # convert to sir Kelvin
    tempd.convert_units('K') # convert to kelvin
    
    press0 = press.interpolate([('longitude',[latlon[1]]),
                                ('latitude',[latlon[0]])],
                               iris.analysis.Linear())
    temp0  = temps.interpolate([('longitude',[latlon[1]]),
                                ('latitude',[latlon[0]])],
                               iris.analysis.Linear())
    tempd0 = tempd.interpolate([('longitude',[latlon[1]]),
                                ('latitude',[latlon[0]])],
                               iris.analysis.Linear())
    
    f160ax = plotting.ax_skewt()
    f160ax.semilogy(np.squeeze(temp0.data),np.squeeze(press0.data), 'k',linewidth=2, label='T')
    f160ax.semilogy(np.squeeze(tempd0.data),np.squeeze(press0.data), 'teal',linewidth=2, label='Td')
    
    
    ## plot a bunch of nearby profiles
    if nearby>0:
        # find index nearest to lat/lon
        lati,loni = utils.lat_lon_index(latlon[0],latlon[1],lats.points,lons.points)
        
        # plot closest, then a range around the closest
        f160ax.semilogy(np.squeeze(temps[:,lati,loni].data), 
                        np.squeeze(press[:,lati,loni].data),
                        alpha=alpha, color='k',
                        label='nearby')
        f160ax.semilogy(np.squeeze(tempd[:,lati,loni].data), 
                        np.squeeze(press[:,lati,loni].data),
                        alpha=alpha, color='teal')
        for i in range(1,nearby+1,1):    
            for latii in [lati-i,lati, lati+i]:
                for lonii in [loni-i, loni, loni+i]:
                    if latii==lati and lonii==loni:
                        continue # don't do middle again
                    f160ax.semilogy(np.squeeze(temps[:,latii,lonii].data), 
                                    np.squeeze(press[:,latii,lonii].data),
                                    alpha=alpha, color='k',
                                    label='nearby')
                    f160ax.semilogy(np.squeeze(tempd[:,latii,lonii].data), 
                                    np.squeeze(press[:,latii,lonii].data),
                                    alpha=alpha, color='teal',
                                    label='nearby')
    # update y ticks (they get scale formatted by semilogy)
    f160ax.yaxis.set_major_formatter(tick.ScalarFormatter())
        
    
def f160_hour(dtime=datetime(2016,1,6,7), latlon = [-32.87,116.1]):
    '''
    Look at F160 plots over time for a particular location
    INPUTS: hour of interest, latlon, and label to describe latlon (for plotting folder)
    '''
    # Use datetime and latlon to determine what data to read
    extentname='sirivan'
    if dtime < datetime(2017,1,1):
        extentname='waroona'
    extent = plotting._extents_[extentname]
    # also for plotname
    latlon_stamp="%.3fS_%.3fE"%(-latlon[0],latlon[1])
    pnames = 'figures/%s/skewt/fig_%s_%s.png'
    
    # read pressure and temperature cubes
    _,_,th1,_= fio.read_waroona(dtime, extent=extent, add_dewpoint=True)#, add_theta=True)
    p,T,Td = th1.extract(['air_pressure','air_temperature','dewpoint_temperature'])
    
    ffdtimes = utils.dates_from_iris(p)
    
    for i in range(len(ffdtimes)):
        dstamp = ffdtimes[i].strftime("%Y%m%d%H%M")
        ptitle="SkewT$_{ACCESS}$   (%s) %s"%(latlon_stamp,ffdtimes[i].strftime("%Y %b %d %H:%M (UTC)"))
        pname=pnames%(extentname,latlon_stamp,dstamp)
        
        plt.close()
        f160(p[i],T[i],Td[i], latlon)
        plt.title(ptitle)
        plt.savefig(pname)
        print("INFO: Saved figure ",pname)
    
if __name__ == '__main__':
    
    print("INFO: testing cloud_outline.py")
    #emberstorm_clouds(datetime(2016,1,5,15))
    
    for dtime in [ datetime(2016,1,6,7) + timedelta(hours=x) for x in range(2) ]:
        #emberstorm_clouds(dtime)
        f160_hour(dtime)#,old=True)

