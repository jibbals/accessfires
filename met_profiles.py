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

import matplotlib
matplotlib.use('Agg')

from metpy.units import units
#distance = np.arange(1, 5) * units.meters # easy way to add units (Also units.metre works!)
#g = 9.81 * units.meter / (units.second * units.second)
from metpy.plots import SkewT,Hodograph

# plotting stuff
import matplotlib.pyplot as plt
#import matplotlib.colors as col
#import matplotlib.ticker as tick
#import matplotlib.gridspec as gridspec

import numpy as np
from datetime import datetime,timedelta

# local modules
from utilities import plotting, utils, fio, constants
from pyrocb import __PCB_occurrences__


###
## GLOBALS
###
_sn_ = 'met_profiles'

def hodograph(u,v,latlon,average=0,ax=None,ztop=18000,axlims=[-50,50],markers=[1000,2000,3000,4000],**hodoargs):
    """
    hodograph plot
    args:
        u: iris cube [lev,lat,lon] east-west wind in m/s
        v: south-north wind
        latlon: [lat,lon] to be profiled
        average (km): interpolate? or use area average? 
        ax: plot axis for hodograph creation
        ztop: only look up to this height (m)
        further arguments will be sent to Hodograph plot function 
    """
    
    height = utils.height_from_iris(u)
    zinds=height > -999
    if ztop is not None:
        zinds=height<ztop
    height=height/1000.0 * units("kilometre") # add units for metpy
    u0 = utils.profile_interpolation(u,latlon,average=average).data * units('m/s')
    v0 = utils.profile_interpolation(v,latlon,average=average).data * units('m/s')
    
    ## Default hodograph args:
    #if 'cmap' not in hodoargs:
    if ax is None:
        hodo=Hodograph()
    else:
        hodo=Hodograph(ax)
    # Default colormap to use to show height
    if 'cmap' not in hodoargs:
        hodoargs['cmap']='nipy_spectral'
    
    cs=hodo.plot_colormapped(u0[zinds],v0[zinds],height[zinds],**hodoargs)
    # add markers at specified heights
    if markers is not None:
        heights = utils.height_from_iris(u) # in metres
        # indices for nearest points
        markinds=np.ones(len(markers),dtype=np.int)
        for mi, mheight in enumerate(markers):
            markinds[mi] = np.abs(heights - mheight).argmin()
        # plot markers to show the height clearly
        hodo.plot(u0[markinds],v0[markinds],
                  marker='x',color='k', linestyle='None')
    hodo.add_grid(increment=20)
    hodo.ax.set_xlim(axlims[0],axlims[1])
    hodo.ax.set_ylim(axlims[0],axlims[1])
    return cs
    
def f160(press,Temp,Tempd, latlon,
         fig=None,
         press_rho=None,uwind=None,vwind=None, 
         nearby=2, alpha=0.3,
         **skewtargs):
    '''
    show skewt logp plot of temperature profile at latlon
    input cubes: pressures, temperatures are [z,lat,lon]
    profile for cubes interpolated to latlon point will be shown with linewidth of 2
    nearest values to lat lon within <nearby> indices will be plotted at low alpha
    
    Wind barbs will be shown if p_ro,u,v are set
    '''
    tcolor='k'
    tdcolor='m'
    
    ## first interpolate pressure and temperature to latlon
    lons=press.coord('longitude')
    lats=press.coord('latitude')
    
    # Extract data array from cube, convert to desired units
    fromunits=[units(str(press.units)), units(str(Temp.units)), units(str(Tempd.units))]
    tounits=[units.mbar, units.degC, units.degC]
    pTTd = []
    for cube,funit,tunit in zip([press,Temp,Tempd],fromunits,tounits):
        # interpolate
        darr0 = utils.profile_interpolation(cube,latlon).data
        # unmask if necessary
        darr1 = utils.unmask(darr0)
        # conversion
        darr2 = (darr1*funit).to(tunit)
        pTTd.append(darr2)
    p,T,Td = pTTd
    
    # default args to SkewT method of library
    if 'fig' not in skewtargs:
        skewtargs['fig'] = plt.figure(figsize=(9,9))
    if 'rotation' not in skewtargs:
        skewtargs['rotation']=30
    
    skew = SkewT(**skewtargs)
    skew.plot(p,T,tcolor, linewidth=2)
    skew.plot(p,Td,tdcolor, linewidth=2)
    
    if np.any(Td > T):
        print("WARNING: Td > T at some point!?")
        print("     Td:",Td)
        print("      T:",T)
    
    ## Add wind profile if desired
    if uwind is not None and vwind is not None and press_rho is not None:
        # interpolate to desired lat/lon
        u0 = utils.unmask(utils.profile_interpolation(uwind,latlon).data)
        v0 = utils.unmask(utils.profile_interpolation(vwind,latlon).data)
        p_rho0 = utils.unmask(utils.profile_interpolation(press_rho,latlon).data)
        
        u = u0 * units('m/s')#units(str(uwind.units))
        v = v0 * units('m/s')#units(str(vwind.units))
        u = u.to(units.knots)
        v = v.to(units.knots)
        pro = (p_rho0 * units(str(press_rho.units))).to(units.mbar)
        #print("DEBUG: f160 interp2", u.shape, v.shape, pro.shape)
        #print("DEBUG: f160 interp2", u,v,pro)
        nicer_z=np.union1d(np.union1d(np.arange(0,41,5), np.arange(43,81,3)), np.arange(81,140,1))
        if latlon[1]>120: # sirivan has more low levels
            nicer_z=np.union1d(np.union1d(np.arange(0,41,8), np.arange(43,81,5)), np.arange(81,140,2))
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
    skew.ax.set_ylim(1000,300)
    skew.ax.set_xlim(-30,60)
    # Add the relevant special lines
    skew.plot_dry_adiabats()
    skew.plot_moist_adiabats()
    skew.plot_mixing_lines()
    return skew

def model_metpy_hour(dtime=datetime(2016,1,6,7), 
                     latlon=constants.latlons['pyrocb_waroona1'],
                     latlon_stamp=None,
                     model_version='waroona_run1',
                     nearby=2,
                     HSkip=None,
                     ztop=15000,
                     ):
    '''
    Look at F160 plots over time for a particular location
    INPUTS: hour of interest, latlon, and label to describe latlon (for plotting folder)
        nearby: how many gridpoints around the desired latlon to plot (can be 0)
    '''
    
    # Use datetime and latlon to determine what data to read
    extentname=model_version.split('_')[0]
    extent = [latlon[1]-0.05, latlon[1]+0.05, latlon[0]-0.05, latlon[0]+0.05]
    
    if latlon_stamp is None:
        latlon_stamp="%.3fS_%.3fE"%(-latlon[0],latlon[1])
    
    # Show marker to make sense of profile latlon
    f,ax = plotting.map_tiff_qgis(extentname+".tiff",)
    plotting.map_add_locations_extent(extentname,nice=True)
    plotting.map_add_nice_text(ax,[latlon],texts=[latlon_stamp],markercolors=['red'])
    fio.save_fig(model_version,_sn_,"map",plt,subdir=latlon_stamp)
    
    # read pressure and temperature cubes
    cubes = fio.read_model_run(model_version, fdtime=dtime, extent=extent,
                               add_winds = True,
                               add_dewpoint = True,
                               HSkip=HSkip)
    p, t, td, u, v = cubes.extract(['air_pressure','air_temperature',
                                    'dewpoint_temperature','u','v'])
    
    p_rho = p
    if model_version=='waroona_run1':
        p_rho, = cubes.extract('air_pressure_rho')
    
    ffdtimes = utils.dates_from_iris(p)
    offset_hours=fio.run_info[model_version]['UTC_offset']
    offset=timedelta(hours=offset_hours)
    
    for i in range(len(ffdtimes)):
        utc=ffdtimes[i]
        lt =utc+offset
        
        ## Plot side by side: FAILS for no reason... axes get messed up
        
        ## f160 plot
        suptitle="%s (%s) %s"%(model_version,latlon_stamp,lt.strftime("%H:%M (LT)"))
        
        
        #print(ffdtimes[i].strftime("DEBUG: %d %H:%M"))
        f160(p[i],t[i],td[i], latlon,
             #fig=fig,subplot=ax_f160,
             press_rho=p_rho[i], uwind=u[i], vwind=v[i],
             )
        
        # Save figure
        plt.title(suptitle)
        fio.save_fig(model_version,_sn_,utc,plt,subdir=latlon_stamp+'/f160')
        
        ## Plot HODOR GRAPHIC
        markers=[1000,2000,4000]
        cs=hodograph(u[i,],v[i,],latlon,ztop=ztop,axlims=[-40,40],markers=markers)
        plt.title(suptitle)
        
        cbar_ax = plt.gcf().add_axes([0.79, 0.3, 0.03, 0.4])# X Y Width Height
        cb = plt.colorbar(cs, cax=cbar_ax, orientation='vertical',
                          format=matplotlib.ticker.ScalarFormatter(), 
                          pad=0)
        cb.set_label('Altitude [km]')
        # add markers to colorbar, axes are 0-1, 0-1
        cb.ax.plot([0.5]*len(markers), np.array(markers)/ztop, 
                   marker='x', color='k',linestyle='None')
        
        # Save figure
        fio.save_fig(model_version,_sn_,utc,plt,subdir=latlon_stamp+'/hodo')
        

if __name__ == '__main__':
    
    w_run = 'waroona_run3'
    w_upwind = []
    w_pcb = [-32.9,116.05] # latlon
    
    
    if True: # lets look at sirivan
        si_runs = ['sirivan_run6_hr',]#'sirivan_run6_hr', ]#'sirivan_run7_hr', 'sirivan_run7']
        for si_run in si_runs:
            
            si_hours = fio.run_info[si_run]['filedates'][8:]
            si_pcb = __PCB_occurrences__[si_run]['latlon'][-1]
            si_upwind = constants.latlons['fire_sirivan_upwind']
            si_mid = -32, 149.8 # sirivan middle of burn area
            
            for hour in si_hours:
                for latlon,stamp in zip([si_upwind,si_mid,si_pcb],
                                        ['upwind', 'mid_fire', 'pcb']):
                    model_metpy_hour(dtime=hour,
                                     latlon=latlon,
                                     latlon_stamp=stamp,
                                     model_version=si_run, 
                                     nearby=0,
                                     HSkip=None)
    
