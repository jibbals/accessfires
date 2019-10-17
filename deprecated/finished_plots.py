# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 10:11:25 2019
    
    old plots that I don't know if I need any more
    
@author: jgreensl
"""

import matplotlib
# don't plot on screen, send straight to file
# this is for lack of NCI display
matplotlib.use('Agg')

# plotting stuff
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np

# local modules
from utilities import plotting, utils, fio


####################################
############ METHODS ###############
####################################
    
    
def jeff_transect(pname="figures/jeff_transect_sample.png"):
    
    nlev=70 # Save on ram just look at lower levels
    tstep=0
    
    dat = fio.read_pcfile('data/umnsaa_pc2016010515.nc')
    lat  = dat['latitude'   ][:]
    lon  = dat['longitude'  ][:]
    w   = dat['upward_air_velocity'][tstep,:nlev,:,:]
    qc  = dat['mass_fraction_of_cloud_liquid_water_in_air'][tstep,:nlev,:,:] + dat['mass_fraction_of_cloud_ice_in_air'][tstep,:nlev,:,:]
    # Also some stuff based on calculated data (in fio.py)
    theta = dat['theta'][tstep,:nlev,:,:]
    zth = dat['z_theta'][tstep,:nlev,:,:]
    s   = dat['wind_speed'][tstep,:nlev,:,:]
    # Save some ram:
    del dat
    nz,ny,nx = qc.shape
    # READ TOPOG DATA FROM PA
    topog, latt, lont = fio.read_topog('data/umnsaa_pa2016010515.nc')
    
    # cross section lat,lon start and finish
    transects = plotting._transects_
    start,end = transects['waroona1']
    npoints = 50
    
    # Pull out cross section of topography and height
    slicetopog = utils.cross_section(topog,start,end, latt,lont,npoints=npoints)
    slicez = utils.cross_section(zth,start,end,lat,lon,npoints=npoints)
    xticks,xlabels = utils.cross_section_ticks_labels(start,end)
    xaxis=np.linspace(0,1,npoints)
    slicex=np.tile(xaxis,(nz,1))
    
    ## Set up figure of 2x2
    f,axes = plt.subplots(2,2, sharex=True, sharey=True)
    
    # Potential temperature
    slicetheta = utils.cross_section(theta,start,end,lat,lon,npoints=npoints)
    thetalevels = np.arange(280,320,2) # contour lines to plot
    thetacontours = thetalevels
    thetacmap = plt.cm.get_cmap('YlOrRd')
    
    # Vertical velocity
    slicew = utils.cross_section(w,start,end,lat,lon,npoints=npoints)
    wlevels = 2.0**np.arange(-2,6)
    wlevels = np.union1d(np.union1d(wlevels,-wlevels),np.array([0]))
    wcontours = np.array([0])
    cmapw = plt.cm.get_cmap('PiYG')
    cmapw.set_over('k')
    wnorm = col.SymLogNorm(0.25)
    
    # wind speed
    s[np.isnan(s)] = -5000 # There is one row or column of s that is np.NaN, one of the edges I think
    slices = utils.cross_section(s,start,end,lat,lon,npoints=npoints)
    slevels = np.arange(0,22,2)
    cmaps = plt.cm.get_cmap('YlGnBu')
    
    # clouds by water+ice content
    sliceqc = utils.cross_section(qc,start,end,lat,lon,npoints=npoints)
    qccontours = np.arange(0.0,2.25,0.25)
    cmapqc = plt.cm.get_cmap('YlGnBu')
    
    ztop=4000
    for ax, slicedata, slicelevels, cmap, slicecontours, title,norm, cbarform in [
        [axes[0,0], slicetheta,thetalevels,thetacmap,thetacontours, "T$_{\\theta}$ (K)", None, None],
        [axes[0,1], slicew, wlevels, cmapw, wcontours, "Vertical motion (m/s)",wnorm, tick.ScalarFormatter()],
        [axes[1,0], slices, slevels, cmaps, slevels, "Wind (m/s)", None, None],
        [axes[1,1], sliceqc, qccontours, cmapqc, qccontours, "Water+ice (kg/kg air)",None,None]
        ]:
        plt.sca(ax)
        # Note that contourf can work with non-plaid coordinate grids provided both are 2-d
        # Contour inputs: xaxis, yaxis, data, colour gradient 
        plt.contourf(slicex,slicez,slicedata,slicelevels,cmap=cmap,norm=norm)
        plt.colorbar(format=cbarform)
        # Add contour lines
        plt.contour(slicex,slicez,slicedata,slicecontours,colors='k')            
        # make sure land is obvious
        plt.fill_between(xaxis,slicetopog,interpolate=True,facecolor='black')
        plt.xticks([])
        if ztop != None:
            plt.ylim(0,ztop)
        plt.xlabel('')
        plt.title(title)
    
    for ax in [axes[0,0],axes[1,0]]:
        plt.sca(ax)
        plt.ylabel('Height (m)')
    
    for ax in [axes[1,0],axes[1,1]]:
        plt.sca(ax)
        #plt.xticks(xticks,xlabels,rotation=-15)
        plt.xlabel("Transect")
        #plt.xticks(xticks,xlabels, rotation=25)
        plt.xticks([])
    
    
    plt.savefig(pname)
    print("FIGURE SAVED: ",pname)
    plt.close()
    
if __name__ == '__main__':
    
    jeff_transect()
