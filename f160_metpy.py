# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 12:28:05 2019
    F160 plot using MetPy
@author: jgreensl
"""

## metpy imports
## conda install -c conda-forge metpy
## 
## NO WAY! Also contained in NCI
## module use /g/data3/hh5/public/modules
## module load conda/analysis3

import metpy
from metpy.units import units
#distance = np.arange(1, 5) * units.meters # easy way to add units (Also units.metre works!)
#g = 9.81 * units.meter / (units.second * units.second)
from metpy.plots import SkewT

import matplotlib
#matplotlib.use('Agg',warn=False)

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

def f160(press,Temp,Tempd, latlon, 
         press_ro=None,uwind=None,vwind=None, 
         nearby=2, alpha=0.3):
    '''
    show skewt logp plot of temperature profile at latlon
    input cubes: pressures, temperatures are [z,lat,lon]
    profile for cubes interpolated to latlon point will be shown with linewidth of 2
    nearest values to lat lon within <nearby> indices will be plotted at low alpha
    
    Wind barbs will be shown if p_ro,u,v are set
    '''
    tcolor='r'
    tdcolor='g'
    
    ## first interpolate pressure and temperature to latlon
    lons=press.coord('longitude')
    lats=press.coord('latitude')
    #press.convert_units('hPa') # convert to hPa
    #Temp.convert_units('C') # convert to sir Kelvin
    #Tempd.convert_units('C') # convert to kelvin
    
    press0 = press.interpolate([('longitude',[latlon[1]]),
                                ('latitude',[latlon[0]])],
                               iris.analysis.Linear())
    temp0  = Temp.interpolate([('longitude',[latlon[1]]),
                            ('latitude',[latlon[0]])],
                           iris.analysis.Linear())
    tempd0 = Tempd.interpolate([('longitude',[latlon[1]]),
                             ('latitude',[latlon[0]])],
                            iris.analysis.Linear())
    
    # Plot T, and Td
    # pull out data array (units don't work with masked arrays)
    p = np.squeeze(press0.data.data) * units(str(press.units))
    T = np.squeeze(temp0.data.data) * units(str(Temp.units))
    T = T.to(units.degC)
    Td = np.squeeze(tempd0.data.data) * units(str(Tempd.units))
    Td = Td.to(units.degC)
    
    fig = plt.figure(figsize=(9,9))
    skew = SkewT(fig,rotation=45)
    skew.plot(p,T,tcolor, linewidth=2)
    skew.plot(p,Td,tdcolor, linewidth=2)
    
    
    ## Add wind profile if desired
    if uwind is not None and vwind is not None and press_ro is not None:
        # interpolate to desired lat/lon
        u0 = uwind.interpolate([('longitude',[latlon[1]]),
                            ('latitude',[latlon[0]])],
                           iris.analysis.Linear())
        v0 = vwind.interpolate([('longitude',[latlon[1]]),
                            ('latitude',[latlon[0]])],
                           iris.analysis.Linear())
        p_ro0 = press_ro.interpolate([('longitude',[latlon[1]]),
                            ('latitude',[latlon[0]])],
                           iris.analysis.Linear())
        u   = np.squeeze(u0.data.data) * units('m/s')#units(str(uwind.units))
        v   = np.squeeze(v0.data.data) * units('m/s')#units(str(vwind.units))
        u   = u.to(units.knots)
        v   = v.to(units.knots)
        pro = np.squeeze(p_ro0.data.data) * units(str(press_ro.units))
        nicer_z=np.union1d(np.union1d(np.arange(0,41,5), np.arange(43,81,3)), np.arange(81,140,1))
        #skip=(slice(None,None,None),nicer_z)
        skew.plot_barbs(pro[nicer_z],u[nicer_z],v[nicer_z])
        
    ## plot a bunch of nearby profiles
    if nearby>0:
        # find index nearest to lat/lon
        lati,loni = utils.lat_lon_index(latlon[0],latlon[1],lats.points,lons.points)
        
        for i in range(1,nearby+1,1):    
            for latii in [lati-i,lati, lati+i]:
                for lonii in [loni-i, loni, loni+i]:
                    # plot a range around the closest
                    p = np.squeeze(press[:,latii,lonii].data.data) * units(str(press.units))
                    T = np.squeeze(Temp[:,latii,lonii].data.data) * units(str(Temp.units))
                    T = T.to(units.degC)
                    Td = np.squeeze(Tempd[:,latii,lonii].data.data) * units(str(Tempd.units))
                    Td = Td.to(units.degC)
                    skew.plot(p,T,tcolor, linewidth=1, alpha=alpha)
                    skew.plot(p,Td,tdcolor, linewidth=1, alpha=alpha)

    # set limits
    skew.ax.set_ylim(1000,100)
    skew.ax.set_xlim(-30,60)
    # Add the relevant special lines
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()
    return skew

def f160_hour(dtime=datetime(2016,1,6,7), latlon=plotting._latlons_['pyrocb_waroona1'],
              subfolder='metpy_pyrocb_waroona1',old=False):
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
        p,t,td = cubes.extract(['air_pressure','air_temperature','dewpoint_temperature'])
        pro, u, v = cubes.extract(['air_pressure','u','v'])
    else:
        _,ro1,th1,_= fio.read_waroona(dtime, extent=extent, add_dewpoint=True, add_winds=True)#, add_theta=True)
        p,t,td = th1.extract(['air_pressure','air_temperature','dewpoint_temperature'])
        pro, u, v  = ro1.extract(['air_pressure','u','v'])
    
    ffdtimes = utils.dates_from_iris(p)
    
    
    
    for i in range(len(ffdtimes)):
        # Plot name and title
        dstamp = ffdtimes[i].strftime("%Y%m%d%H%M")
        ptitle="SkewT$_{ACCESS}$   (%s) %s"%(latlon_stamp,ffdtimes[i].strftime("%Y %b %d %H:%M (UTC)"))
        pname=pnames%(extentname,subfolder,latlon_stamp,dstamp)
        
        # create plot
        
        f160(p[i],t[i],td[i], latlon,
             press_ro=pro[i], uwind=u[i], vwind=v[i])
        plt.title(ptitle)
        
        # save plot
        fio.save_fig(pname,plt)
        plt.close()

if __name__ == '__main__':
    
    topleft=[-32.75, 115.8] # random point away from the fire influence
    for dtime in [ datetime(2016,1,6,3) + timedelta(hours=x) for x in range(6) ]:
        f160_hour(dtime,subfolder='metpy_pyrocb_waroona1',old=False)
        f160_hour(dtime,subfolder='metpy_pyrocb_waroona1_old',old=True)
        f160_hour(dtime,latlon=topleft,subfolder='metpy_topleft_waroona',old=False)
        f160_hour(dtime,latlon=topleft,subfolder='metpy_topleft_waroona_old',old=True)

