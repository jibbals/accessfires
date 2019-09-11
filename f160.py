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


def f160(press,temps, tempd, latlon, p_ro=None,u=None,v=None, nearby=2, alpha=0.4):
    '''
    show skewt logp plot of temperature profile at latlon
    Does model give us dewpoint or equiv to also plot here??
    input cubes: pressures, temperatures are [z,lat,lon]
    profile for cubes interpolated to latlon point will be shown with linewidth of 2
    nearest values to lat lon within <nearby> indices will be plotted at low alpha
    
    Wind barbs will be shown if p_ro,u,v are set
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
    
    ## Add wind profile if desired
    if u is not None and v is not None and p_ro is not None:
        # interpolate to desired lat/lon
        p_ro.convert_units('hPa') # convert to hPa
        u0 = u.interpolate([('longitude',[latlon[1]]),
                            ('latitude',[latlon[0]])],
                           iris.analysis.Linear())
        v0 = v.interpolate([('longitude',[latlon[1]]),
                            ('latitude',[latlon[0]])],
                           iris.analysis.Linear())
        p_ro0 = p_ro.interpolate([('longitude',[latlon[1]]),
                            ('latitude',[latlon[0]])],
                           iris.analysis.Linear())
        n_z, _, _ = p_ro0.shape
        ## show profile of horizontal winds on the right
        barbax = plt.twiny(f160ax)
        plt.sca(barbax)
        # plt.barbs([X,Y],U,V) # X,Y are barb locations, U,V are directions
        # all X locations are 1, Y locations should match pressure levels
        # 2D(surface) plot but just making a single line of barbs
        X=np.ones(press[:,0,0].shape);Y=np.squeeze(p_ro0.data)
        U=np.squeeze(u0.data); V=np.squeeze(v0.data)
        mX,mY = np.meshgrid(X,Y)
        mV = np.repeat(V[np.newaxis,:],n_z,0) # repeat along unused x direction
        mV.mask = True
        mV.mask[0,:]=False
        mU = np.repeat(U[np.newaxis,:],n_z,0)
        mU.mask = True
        mU.mask[0,:]=False
        # make it not so crowded at low altitudes
        nicer_z=np.union1d(np.union1d(np.arange(0,41,5), np.arange(43,81,3)), np.arange(81,140,1))
        skip=(slice(None,None,None),nicer_z)
        plt.barbs(mX[skip],mY.transpose()[skip],mU[skip],mV[skip])
        # remove additional new x axis
        plt.xticks([],[])
        plt.xlim([0,1.05])
        plt.sca(f160ax)
    
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
    plt.yticks(np.linspace(1000,100,10))
    #f160ax.xaxis.set_major_locator(tick.MultipleLocator(10))
    
def f160_hour(dtime=datetime(2016,1,6,7), latlon=plotting._latlons_['pyrocb_waroona1'],
              subfolder='pyrocb_waroona1',old=False):
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
    pnames = 'figures/%s/skewt/%s/fig_%s_%s.png'
    
    # read pressure and temperature cubes
    if old:
        # I don't know how to calc ro levels, but should be similar I think
        cubes = fio.read_waroona_pcfile(dtime,extent=extent, add_winds=True, add_dewpoint=True)
        p,T,Td = cubes.extract(['air_pressure','air_temperature','dewpoint_temperature'])
        pro, u, v = cubes.extract(['air_pressure','u','v'])
    else:
        _,ro1,th1,_= fio.read_waroona(dtime, extent=extent, add_dewpoint=True, add_winds=True)#, add_theta=True)
        p,T,Td  = th1.extract(['air_pressure','air_temperature','dewpoint_temperature'])
        pro,u,v = ro1.extract(['air_pressure','u','v'])
    
    ffdtimes = utils.dates_from_iris(p)
    
    for i in range(len(ffdtimes)):
        # Plot name and title
        dstamp = ffdtimes[i].strftime("%Y%m%d%H%M")
        ptitle="SkewT$_{ACCESS}$   (%s) %s"%(latlon_stamp,ffdtimes[i].strftime("%Y %b %d %H:%M (UTC)"))
        pname=pnames%(extentname,subfolder,latlon_stamp,dstamp)
        # create plot
        f160(p[i],T[i],Td[i], latlon,p_ro=pro[i], u=u[i], v=v[i])
        plt.title(ptitle)
        # save plot
        fio.save_fig(pname,plt)
        plt.close()
    
if __name__ == '__main__':
    
    print("INFO: testing cloud_outline.py")
    #emberstorm_clouds(datetime(2016,1,5,15))
    
    for dtime in [ datetime(2016,1,6,3) + timedelta(hours=x) for x in range(8) ]:
        f160_hour(dtime,subfolder='pyrocb_waroona1',old=False)
        f160_hour(dtime,subfolder='pyrocb_waroona1_old',old=True)

