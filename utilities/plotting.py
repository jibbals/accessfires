#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 14:06:49 2019

  In this script is repeatable generic plotting stuff like adding a scale to a cartopy map

@author: jesse
"""

# Plotting stuff
import numpy as np
import matplotlib
from matplotlib import patheffects
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib.projections import register_projection
from matplotlib.collections import LineCollection

# Mapping stuff
import cartopy
import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt

# local stuff
from .context import utils, constants, skewt, thermo

###
### GLOBALS
###

## Standard colour maps for use between figures
_cmaps_ = {}
_cmaps_['verticalvelocity'] = 'PiYG_r'  # jeff likes this for vert velocity
_cmaps_['windspeed']        = 'YlGnBu' # jeff likes this for wind speeds
_cmaps_['qc']               = 'BuPu' # white through to purple for water+ice content in air
_cmaps_['topog']            = 'terrain'
# Extents: EWSN
_extents_               = {}
_latlons_               = {}

# Waroona locations
_extents_['waroona']    = [115.775,116.2, -33.05,-32.7] # local
_extents_['waroonas']   = [112,120,-34.5,-31] # synoptic
_latlons_['waroona']    = -32.84, 115.93  # suburb centre: -32.8430, 115.8526
_latlons_['yarloop']    = -32.96, 115.90  # suburb centre: -32.9534, 115.9124
_latlons_['perth']      = -31.9505, 115.8605
_latlons_['fire_waroona'] = -32.89, 116.17
# two PyroCB
_latlons_['pyrocb_waroona1'] = -32.87,116.1 # ~4pm first day
_latlons_['pyrocb_waroona2'] = 0,0 # 1100-1400 second day

# Sir Ivan locations
_extents_['sirivan']    = [149.2, 150.4, -32.4, -31.6]
_extents_['sirivans']   = [145,154, -34, -30]
_latlons_['dunedoo']    = -31.99, 149.53
_latlons_['uarbry']      = -32.047280, 149.71
_latlons_['cassillis']      = -32.01, 150.0
_latlons_['fire_sirivan'] = -32.45, 149.51
# one pyrocb
_latlons_['pyrocb_sirivan'] = 0,0 # 0530UTC=XXXX local time
# The pycb over Sir Ivan was around 0530 UTC on Sunday 12 Feb (or a bit after).
# The exact location is a bit tricky, because the pycb would be downstream of 
# the actual fire front, but around the location of Uarbry 
# (which was burned over) is probably a good place to start. 
# You'll probably have to move the cross section around a bit to see what 
# lat/longs get the most interesting picture.


_transects_             = {} 
__x0__,__x1__ = 115.8, 116.19

_transects_['waroona1'] = [-32.79   , __x0__], [-32.92   , __x1__]
_transects_['waroona2'] = [-32.82   , __x0__], [-32.93   , __x1__]
_transects_['waroona3'] = [-32.86   , __x0__], [-32.88   , __x1__]
# again but start more southerly and end more northerly
_transects_['waroona4'] = [-32.92   , __x0__], [-32.82   , __x1__]
_transects_['waroona5'] = [-32.96   , __x0__], [-32.85   , __x1__] 
_transects_['waroona6'] = [-32.87   , __x0__], [-32.89   , __x1__]
_transects_['pyrocb_waroona'] = [-32.89 , 115.9 ], [-32.88   , __x1__]
_transects_['pyrocbx1_waroona'] = [-32.95, 116.09], [-32.78   , 116.16]
_transects_['pyrocbx2_waroona'] = [-32.95, 116.15], [-32.78   , 116.09]

# looking at sir ivan
__si0__, __si1__ = 149.4, 149.9
_transects_['sirivan1'] = [-32.  , __si0__], [-32.3 , __si1__]
_transects_['sirivan2'] = [-32.05, __si0__], [-32.35, __si1__] # like sirivan1 but lower
_transects_['sirivan3'] = [-31.95, __si0__], [-32.25, __si1__]   # '' but higher
_transects_['sirivan4'] = [-32.25, __si0__], [-31.95, __si1__] # low lat to high lat crosssec
_transects_['sirivan5'] = [-32.35, __si0__], [-32.05, __si1__] 
_transects_['sirivan6'] = [-32.3 , __si0__], [-32   , __si1__] 

def init_plots():
  matplotlib.rcParams['font.size'] = 15.0
  matplotlib.rcParams["text.usetex"]      = False     # I forget what this is for, maybe allows latex labels?
  matplotlib.rcParams["legend.numpoints"] = 1         # one point for marker legends
  matplotlib.rcParams["figure.figsize"]   = (9, 7)    # Default figure size
  matplotlib.rcParams["axes.titlesize"]   = 19        # title font size
  matplotlib.rcParams["figure.titlesize"] = 21        # figure suptitle size
  matplotlib.rcParams["axes.labelsize"]   = 15        #
  matplotlib.rcParams["xtick.labelsize"]  = 11        #
  matplotlib.rcParams["ytick.labelsize"]  = 11        #
  matplotlib.rcParams['image.cmap'] = 'plasma'        # Colormap default
  matplotlib.rcParams['axes.formatter.useoffset'] = False # another one I've forgotten the purpose of
  # rcParams["figure.dpi"] 
  #matplotlib.rcParams["figure.dpi"] = 400           # DEFAULT DPI for plot output
  # THIS MESSES UP THE PLOTTING

def ax_skewt(tlims=[240,330],plims=[1050,100], th_bl=None, q_bl=None):
    '''
    Using Jeff's saturationpoint code, create base for log p skew t plot
    EG: 
        ax=ax_skewt(tlims=[250,375],plims=[1050,50])
        ax.semilogy(T[0,:,50,50].data, p[0,:,50,50].data/100., 'k')
        plt.show()
    Inputs: 
        temperature limits, pressure limits for plotting
    Returns:
        axis with stuff already drawn
    Additional: 
        boundary layer theta and mixing ratio can be added to draw SP curve 
        
    
    '''
    
    #import matplotlib as mpl
    #import matplotlib.pyplot as plt
    #from matplotlib.ticker import ScalarFormatter#, MultipleLocator
    #from matplotlib.projections import register_projection
    #import numpy as np
    
    #import skewt
    #import thermo
    
    # Dry adiabats rise from right to left, representing constant potential temperature
    def plot_dryadiabat(ax,theta,p,p0):
        # p should be an array of pressures at which the adiabat is plotted
        T = theta*(p/p0)**thermo.kappa
        ax.semilogy(T, p, 'g-',alpha=0.5)
        
    # mixing ratios rise left to right
    def plot_mixrat(ax,r,p):
        # p should be an array of pressures (hPa) at which the adiabat is plotted
        e = p * r / (thermo.eps + r)
        loge = np.log(e)
        Td = (243.5*loge - 440.8)/(19.48 - loge) + 273.15
        ax.semilogy(Td,p,'g-',alpha=0.5)
    
    def plot_moist_adiabats(ax, t0=None, p=None, **kwargs):
        r"""Plot moist adiabats.
        
        Temperatures in K, pressures in hPa.

        Adds saturated pseudo-adiabats (lines of constant equivalent potential
        temperature) to the plot. The default style of these lines is dashed
        blue lines with an alpha value of 0.5. These can be overridden using
        keyword arguments.

        Parameters
        ----------
        t0 : array_like, optional
            Starting temperature values in Kelvin. If none are given, they will be
            generated using the current temperature range at the bottom of
            the plot.
        p : array_like, optional
            Pressure values to be included in the moist adiabats. If not
            specified, they will be linearly distributed across the current
            plotted pressure range.
        kwargs
            Other keyword arguments to pass to :class:`matplotlib.collections.LineCollection`

        Returns
        -------
        matplotlib.collections.LineCollection
            instance created

        See Also
        --------
        :func:`~metpy.calc.thermo.moist_lapse`
        :meth:`plot_dry_adiabats`
        :class:`matplotlib.collections.LineCollection`

        """
        # Determine set of starting temps if necessary
        if t0 is None:
            xmin, xmax = ax.get_xlim()
            t0 = np.arange(240, xmax, 10)
        print("DEBUG: t0",t0.shape, t0)
        # Get pressure levels based on ylims if necessary
        if p is None:
            p = np.linspace(1050,100, 96)
        print("DEBUG: p",p.shape, p)
        # Assemble into data for plotting
        t = thermo.moist_lapse(p*100, t0[:, np.newaxis], 100000) # in Kelvin
        linedata = [np.vstack((ti, p)).T for ti in t]
        print("DEBUG: t",t.shape, t)
        # Add to plot
        kwargs.setdefault('colors', 'b')
        kwargs.setdefault('linestyles', 'dashed')
        kwargs.setdefault('alpha', 0.5)
        #return ax.add_collection(LineCollection(linedata, **kwargs))
        ax.add_collection(LineCollection(linedata, **kwargs))
    
    #    def plot_moistadiabat(ax, T, p, p0=1e5):
    #        # moist adiabat, 
    #        #at the saturated mixing ratios (kg/kg) over pressure (Pa)
    #        
    #        #for Tk in np.array(T):
    #        r_sat = thermo.r_sat(T,p0)# kg/kg = f(K,Pa)
    #        th_e = thermo.theta_e(T, r_sat, p0) # K/Pa = f(K, kg/kg, Pa)
    #        print("DEBUG: th_e",th_e.shape, th_e)
    #        print("DEBUG: p",p.shape, p)
    #        print("DEBUG: r_sat",r_sat.shape,r_sat)
    #        e = p * r_sat / (thermo.eps + r_sat) 
    #        loge = np.log( 0.01*np.maximum(e, 1e-100) )
    #        T_lcl = 2840.0 / (3.5*np.log(T[0]) - loge - 4.805) + 55.0
    #        #print("DEBUG:",th_e.shape, p.shape, r_sat.shape, T_lcl.shape)
    #        Te = ( th_e * (100000/p0) ** (.079912*r_sat - thermo.kappa) 
    #              * np.exp((2.54-3376/T_lcl)*r_sat*(1+.81*r_sat)))
    #        print("DEBUG: Te",Te.shape, Te)
    #        ax.semilogy(Te,p/100.,'g:',linewidth=2)
    #        #for Tk in np.array(T):
    #        #    r_sat = thermo.r_sat(Tk,p) # kg/kg = f(K,Pa)
    #        #    
    #        #    th_e = thermo.theta_e(Tk, r_sat, p) # K/Pa = f(K, kg/kg, Pa)
    #        #    e = p * r_sat / (thermo.eps + r_sat) 
    #        #    loge = np.log( 0.01*np.maximum(e, 1e-100) )
    #        #    T_lcl = 2840.0 / (3.5*np.log(Tk) - loge - 4.805) + 55.0
    #        #    #print("DEBUG:",th_e.shape, p.shape, r_sat.shape, T_lcl.shape)
    #        #    Te = ( th_e * (100000/p) ** (.079912*r_sat - thermo.kappa) 
    #        #          * np.exp((2.54-3376/T_lcl)*r_sat*(1+.81*r_sat)))
    #        #    ax.semilogy(Te,p,'g:')
    #        #print("DEBUG: Te",Te)
    #        #print("DEBUG: p",p)
    #        #print("DEBUG: r_sat",r_sat)
    #        #print("DEBUG: T_lcl",T_lcl)
        
    matplotlib.rcParams['font.size'] = 18.0
    matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
    mtfs = 22.0  # fontsize for math text
    
    register_projection(skewt.SkewXAxes)
    
    p0 = 1e5
    
    fig = plt.figure(1)#figsize=(6.5875, 6.2125))
    ax = fig.add_subplot(111, projection='skewx')
    
    yticks = np.ravel(np.outer(10**np.arange(0,3),np.linspace(1,10,10)))
    
    plt.grid(True)
    for thx in range(300,2000,100):
        plot_dryadiabat(ax,thx,yticks,1e3)
    for rx in [1e-4,2e-4,4e-4,1e-3,2e-3,4e-3,1e-2,2e-2,4e-2,0.1,0.2,0.4,1.0]:
        plot_mixrat(ax,rx,yticks)
    plot_moist_adiabats(ax)
    #for Tx in range(270,325,10):
    #    Tarr = np.zeros([20])+Tx
    #    parr = np.logspace(5,4,20)
    #    plot_moistadiabat(ax, Tarr,parr)
    
    
    
    if th_bl is not None and q_bl is not None:
        #th_bl = 303
        #q_bl = 0.005
        
        T_bl = th_bl*(p0/1e5)**(thermo.kappa)
        r_bl = q_bl/(1 - q_bl)
        e_bl = p0 * r_bl / (thermo.eps + r_bl)
        loge_bl = np.log(0.01*e_bl)
        Td_bl = (243.5*loge_bl - 440.8)/(19.48 - loge_bl) + 273.15
    
        gamma = 2.0
        delta = 6.6*1e3  # need the 1e3 to get units in K/(kg/kg)
    
        th_fire = gamma*th_bl
        T_fire = th_fire*(p0/1e5)**(thermo.kappa)
        q_fire = ((gamma-1)/delta)*th_bl + 0.86*q_bl
        r_fire = q_fire/(1 - q_fire)
        e_fire = p0 * r_fire / (thermo.eps + r_fire)
        loge_fire = np.log(0.01*e_fire)
        Td_fire = (243.5*loge_fire - 440.8)/(19.48 - loge_fire) + 273.15
    
        al = np.linspace(0.0,1.0,1000)
    
        th_sp = (1-al)*th_bl + al*th_fire
        q_sp = (1-al)*q_bl + al*q_fire
        r_sp = q_sp/(1 - q_sp)
    
        T0_sp = th_sp*(p0/1e5)**(thermo.Rd/thermo.Cpd)
    
        e = p0 * r_sp / (thermo.eps + r_sp)
        loge = np.log( 0.01*np.maximum(e, 1e-100) )
        T_sp = 2840.0 / (3.5*np.log(T0_sp) - loge - 4.805) + 55.0
        p_sp = 1e5*(T_sp/th_sp)**(thermo.Cpd/thermo.Rd)
    
        # Plot the data using normal plotting functions, in this case using
        # log scaling in Y, as dicatated by the typical meteorological plot
        ax.semilogy(T_sp, p_sp*1e-2, 'r')
        #ax.semilogy([T_fire,T_sp[-1],Td_fire],[p0*1e-2,p_sp[-1]*1e-2,p0*1e-2],'b-')
        ax.semilogy([T_bl,T_sp[0],Td_bl],[p0*1e-2,p_sp[0]*1e-2,p0*1e-2],'b-')
        
        ax.text(0.8,0.9,  r'$\theta_{bl}='  +'{:.0f}'.format(th_bl)   +'\mathrm{K}$',transform=ax.transAxes,fontsize=mtfs)
        ax.text(0.8,0.825,r'$q_{bl}='       +'{:.1f}'.format(q_bl*1e3)+'\mathrm{g/kg}$',transform=ax.transAxes,fontsize=mtfs)
        ax.text(0.8,0.75, r'$\gamma='       +'{:.1f}$'.format(gamma),transform=ax.transAxes,fontsize=mtfs)
        ax.text(0.8,0.675,r'$\delta='       +'{:.1f}$'.format(delta*1e-3),transform=ax.transAxes,fontsize=mtfs)
        ax.text(0.8,0.6,  r'$\theta_{fire}='+'{:.0f}'.format(th_fire)+'\mathrm{K}$',transform=ax.transAxes,fontsize=mtfs)
        ax.text(0.8,0.525,r'$q_{fire}='     +'{:.1f}'.format(q_fire*1e3)+'\mathrm{g/kg}$',transform=ax.transAxes,fontsize=mtfs)
    
    # Disables the log-formatting that comes with semilogy
    ax.yaxis.set_major_formatter(tick.ScalarFormatter())
    ax.set_yticks(yticks)
    ax.set_ylim(plims[0], plims[1])
    
    ax.xaxis.set_major_locator(tick.MultipleLocator(10))
    #ax.set_xlim(Td_bl,400)
    ax.set_xlim(tlims[0],tlims[1])
    
    ax.set_xlabel('$T \, (\mathrm{K})$')
    ax.set_ylabel('$p \, (\mathrm{hPa})$')
    
    fig.set_size_inches(20,10)
    #plt.show()
    #plt.savefig('skewT_th{:.0f}_q{:.0f}_gam{:.0f}_del{:.0f}.png'.format(th_bl,q_bl*1e3,gamma,delta*1e-3),dpi=200)
    return ax 

def map_add_locations_extent(extentname, hide_text=False):
    '''
    
    '''
    if extentname == 'waroona':
        map_add_locations(['waroona','yarloop'], 
                          text=[['Waroona', 'Yarloop'],['','']][hide_text], 
                          textcolor='k')
        # add fire ignition
        map_add_locations(['fire_waroona'],
                          text = [['Fire ignition'],['']][hide_text], 
                          color='r', marker='*', 
                          textcolor='k')
        # add pyroCB
    else:
        map_add_locations(['sirivan','uarbry'], 
                          text=[['Sir Ivan','Uarbry'],['','']][hide_text],
                          dx=[-.02,.05], dy =[-.015,-.03],
                          textcolor='k')
        # add fire ignition
        map_add_locations(['fire_sirivan'],
                          text = [['Fire ignition']['']][hide_text], dx=.05,
                          color='r', marker='*', 
                          textcolor='k')

def transect(data, z, lat, lon, start, end, npoints=100, 
             topog=None, latt=None, lont=None, ztop=4000,
             title="", ax=None, colorbar=True,
             cmap=None, norm=None, cbarform=None,
             contours=None,lines=None):
    '''
    Draw cross section
        data is 3d
        z (3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    # Potential temperature
    slicedata  = utils.cross_section(data,lat,lon,start,end,npoints=npoints)
    
    # Pull out cross section of topography and height
    if latt is None:
        latt=lat
    if lont is None:
        lont=lon
    if topog is not None:
        slicetopog = utils.cross_section(topog,latt,lont,start,end,npoints=npoints)
    
    slicez = utils.cross_section(z,lat,lon,start,end,npoints=npoints)
    #xticks,xlabels = utils.cross_section_ticks_labels(start,end)
    xaxis=np.linspace(0,1,npoints)
    slicex=np.tile(xaxis,(len(z),1))
    
    if ax is not None:
        plt.sca(ax)
    # Note that contourf can work with non-plaid coordinate grids provided both are 2-d
    # Contour inputs: xaxis, yaxis, data, colour gradient 
    if contours is None:
        plt.contourf(slicex,slicez,slicedata,cmap=cmap,norm=norm,extend='max')
    else:
        plt.contourf(slicex,slicez,slicedata,contours,cmap=cmap,norm=norm,extend='max')
    
    if colorbar:
        plt.colorbar(format=cbarform, pad=0.01) # pad is distance from axes
    
    # Add contour lines
    if lines is not None:
        plt.contour(slicex,slicez,slicedata,lines,colors='k')            
    
    # make sure land is obvious
    if topog is not None:
        plt.fill_between(xaxis,slicetopog,interpolate=True,facecolor='black')
    
    if ztop != None:
        plt.ylim(0,ztop)
    
    plt.xticks([])
    plt.xlabel('')
    plt.title(title)

    return slicedata, slicex, slicez

def transect_s(s, z, lat, lon, start, end, npoints=100, 
               topog=None, latt=None, lont=None, ztop=4000,
               title="Wind speed (m/s)", ax=None, colorbar=True,
               cmap='YlGnBu', norm=None, cbarform=None,
               contours=np.arange(0,25,2.5),lines=np.arange(0,25,2.5)):
    '''
    Draw wind speed cross section
        s is 3d wind speed
        z(3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    
    # wind speed
    s[np.isnan(s)] = -5000 # There is one row or column of s that is np.NaN, one of the edges I think
    
    # call transect using some defaults for potential temperature
    return transect(s,z,lat,lon,start,end,npoints=npoints,
                    topog=topog, latt=latt, lont=lont, ztop=ztop,
                    title=title, ax=ax, colorbar=colorbar,
                    cmap=cmap, norm=norm, cbarform=cbarform,
                    contours=contours,lines=lines)

def transect_theta(theta, z, lat, lon, start, end, npoints=100, 
                   topog=None, latt=None, lont=None, ztop=4000,
                   title="$T_{\\theta}$ (K)", ax=None, colorbar=True,
                   cmap='YlOrRd', norm=None, cbarform=None,
                   contours=np.arange(280,320,2),lines=np.arange(280,320,2)):
    '''
    Draw theta cross section
        theta is 3d potential temperature
        z(3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    # call transect using some defaults for potential temperature
    return transect(theta,z,lat,lon,start,end,npoints=npoints,
                    topog=topog, latt=latt, lont=lont, ztop=ztop,
                    title=title, ax=ax, colorbar=colorbar,
                    cmap=cmap, norm=norm, cbarform=cbarform,
                    contours=contours,lines=lines)

def transect_w(w, z, lat, lon, start, end, npoints=100, 
               topog=None, latt=None, lont=None, ztop=4000,
               title="Vertical motion (m/s)", ax=None, colorbar=True,
               cmap=_cmaps_['verticalvelocity'] , norm=col.SymLogNorm(0.25), 
               cbarform=tick.ScalarFormatter(),
               contours=np.union1d(np.union1d(2.0**np.arange(-2,6),-1*(2.0**np.arange(-2,6))),np.array([0])),
               lines=np.array([0])):
    '''
    Draw theta cross section
        w is 3d vertical motion
        z(3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    #if cmap is not None:
    #    if hasattr(cmap,'set_over'):
    #        cmap.set_over('k')
    # call transect using some defaults for vertical velocity w
    return transect(w, z,lat,lon,start,end,npoints=npoints,
                    topog=topog, latt=latt, lont=lont, ztop=ztop,
                    title=title, ax=ax, colorbar=colorbar,
                    cmap=cmap, norm=norm, cbarform=cbarform,
                    contours=contours,lines=lines)

def transect_qc(qc, z, lat, lon, start, end, npoints=100, 
               topog=None, latt=None, lont=None, ztop=4000,
               title="Water and ice (g/kg air)", ax=None, colorbar=True,
               cmap=_cmaps_['qc'] , norm=col.SymLogNorm(0.02), cbarform=tick.ScalarFormatter(),
               contours=np.arange(0.0,0.4,0.01),
               lines=np.array([constants.cloud_threshold])):
    '''
    Draw theta cross section
        qc is 3d vertical motion
        z(3d), lat(1d), lon(1d) is height (m), lats and lons
        start, end are [lat0,lon0], [lat1,lon1]
        contours will be filled colours
        lines will be where to draw black lines
    '''
    
    # call transect using some defaults for vertical velocity w
    return transect(qc, z,lat,lon,start,end,npoints=npoints,
                    topog=topog, latt=latt, lont=lont, ztop=ztop,
                    title=title, ax=ax, colorbar=colorbar,
                    cmap=cmap, norm=norm, cbarform=cbarform,
                    contours=contours,lines=lines)

def transect_ticks_labels(start,end):
  '''
    return xticks and xlabels for a cross section
  '''
  lat1,lon1=start
  lat2,lon2=end
  # Set up a tuple of strings for the labels. Very crude!
  xticks = (0.0,0.5,1.0)
  fmt = '{:.1f}S {:.1f}E'
  xlabels = (fmt.format(-lat1,lon1),fmt.format(-0.5*(lat1+lat2),0.5*(lon1+lon2)),fmt.format(-lat2,lon2))
  return xticks,xlabels


def map_add_locations(namelist, text=None, proj=None,
                      marker='o', color='grey', markersize=None, 
                      textcolor='k',
                      dx=.025,dy=.015):
    '''
    input is list of names to be added to a map using the lat lons in _latlons_
    '''
    for i,name in enumerate(namelist):
        y,x=_latlons_[name]
        # maybe want special text
        if text is not None:
            name = text[i]
        # dx,dy can be scalar or list
        dxi,dyi = dx,dy
        if isinstance(dx, (list,tuple,np.ndarray)):
            dxi,dyi = dx[i], dy[i]

        # Add marker and text
        if proj is None:
            plt.plot(x,y,  color=color, linewidth=0, 
                     marker=marker, markersize=None)
            plt.text(x+dxi, y+dyi, name, color=textcolor,
                     horizontalalignment='right')
        else:
            plt.plot(x,y,  color=color, linewidth=0, 
                     marker=marker, markersize=None, transform=proj)
            plt.text(x+dxi, y+dyi, name, color=textcolor,
                     horizontalalignment='right',
                     transform=proj)
        

    

def map_google(extent,zoom=10,fig=None,subplotxyn=None,draw_gridlines=True,gridlines=None):
    '''
    Draw a map with google image tiles over given EWSN extent
    if this will be part of a larger figure, subplotxyn=[3,3,4] will put it there
    to define where gridlines go, put gridlines = [lats,lons] 
    return figure, axis, projection
    '''
    # Request map from google
    request = cimgt.GoogleTiles()
    gproj=request.crs
    # Use projection to set up plot
    if fig is None:
        fig = plt.figure()
    
    if subplotxyn is None:
        ax = fig.add_subplot(1, 1, 1, projection=gproj)#stamen_terrain.crs)
    else:
        nrows,ncols,n= subplotxyn
        ax = fig.add_subplot(nrows,ncols,n, projection=gproj)

    if draw_gridlines:
        map_draw_gridlines(ax, gridlines=gridlines)

    # Where are we looking
    ax.set_extent(extent)
    # default interpolation ruins location names
    ax.add_image(request, zoom, interpolation='spline36') 
    return fig, ax, gproj

def map_draw_gridlines(ax, linewidth=1, color='black', alpha=0.5, 
                       linestyle='--', draw_labels=True,
                       gridlines=None):
    gl = ax.gridlines(linewidth=linewidth, color=color, alpha=alpha, 
                      linestyle=linestyle, draw_labels=draw_labels)
    if draw_labels:
        gl.xlabels_top = gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
    if gridlines is not None:
        yrange,xrange=gridlines
        gl.xlocator = matplotlib.ticker.FixedLocator(xrange)
        gl.ylocator = matplotlib.ticker.FixedLocator(yrange)

def map_contourf(extent, data, lat,lon, title="",
                 cmap=None, clabel="", clevs=None, norm=None, 
                 cbar=True, cbarform=None):
    '''
    Show topography map matching extents
    '''
    
    cb = plt.contourf(lon,lat,data, levels=clevs, cmap=cmap,norm=norm)
        
    # set x and y limits to match extent
    xlims = extent[0:2] # East to West
    ylims = extent[2:] # South to North
    plt.ylim(ylims); plt.xlim(xlims)
    if cbar:
        cb=plt.colorbar(label=clabel, format=cbarform, pad=0.01)
    plt.title(title)
    ## Turn off the tick values
    plt.xticks([]); plt.yticks([])
    return cb

def map_topography(extent, topog,lat,lon,title="Topography"):
    '''
    Show topography map matching extents
    '''
    # push blue water part of scale a bit lower
    clevs= np.linspace(-150,550,50,endpoint=True)
    cmaptr=plt.cm.get_cmap("terrain")
    return map_contourf(extent, topog,lat,lon,title=title,clevs=clevs,cmap=cmaptr,clabel="m")
    

def utm_from_lon(lon):
    """
    utm_from_lon - UTM zone for a longitude

    Not right for some polar regions (Norway, Svalbard, Antartica)

    :param float lon: longitude
    :return: UTM zone number
    :rtype: int
    """
    return np.floor( ( lon + 180 ) / 6) + 1

def scale_bar(ax, proj, length, location=(0.5, 0.05), linewidth=3,
              units='km', m_per_unit=1000):
    """
    http://stackoverflow.com/a/35705477/1072212
    ax is the axes to draw the scalebar on.
    proj is the projection the axes are in
    location is center of the scalebar in axis coordinates ie. 0.5 is the middle of the plot
    length is the length of the scalebar in km.
    linewidth is the thickness of the scalebar.
    units is the name of the unit
    m_per_unit is the number of meters in a unit
    """
    # find lat/lon center to find best UTM zone
    x0, x1, y0, y1 = ax.get_extent(proj.as_geodetic())
    #x0, x1, y0, y1 = ax.get_extent(proj)
    # Projection in metres
    utm = cartopy.crs.UTM(utm_from_lon((x0+x1)/2))
    # Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(utm)
    # Turn the specified scalebar location into coordinates in metres
    sbcx, sbcy = x0 + (x1 - x0) * location[0], y0 + (y1 - y0) * location[1]
    # Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbcx - length * m_per_unit/2, sbcx + length * m_per_unit/2]
    # buffer for scalebar
    buffer = [patheffects.withStroke(linewidth=5, foreground="w")]
    # Plot the scalebar with buffer
    ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k',
        linewidth=linewidth, path_effects=buffer)
    # buffer for text
    buffer = [patheffects.withStroke(linewidth=3, foreground="w")]
    # Plot the scalebar label
    t0 = ax.text(sbcx, sbcy, str(length) + ' ' + units, transform=utm,
        horizontalalignment='center', verticalalignment='bottom',
        path_effects=buffer, zorder=2)
    left = x0+(x1-x0)*0.05
    # Plot the N arrow
    #t1 = ax.text(left, sbcy, u'\u25B2\nN', transform=utm,
    #    horizontalalignment='center', verticalalignment='bottom',
    #    path_effects=buffer, zorder=2)
    # Plot the scalebar without buffer, in case covered by text buffer
    ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k',
        linewidth=linewidth, zorder=3)
