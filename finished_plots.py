# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 10:11:25 2019
    
    Figures created for AFAC and analysis of sites is saved here
    
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
import matplotlib.patches as mpatches
import numpy as np

# local modules
from utilities import plotting, utils, fio

# maps
import cartopy.crs as ccrs


####################################
############ METHODS ###############
####################################

def clouds_2panel(topog,s,u,v,
                  qc,rh,
                  z,lat,lon,
                  dtime,
                  ffcube=None,
                  extentname='waroona',
                  transect=1, 
                  vectorskip=9,
                  quiverscale=60,
                  ztop=13000,
                  ext='.png',
                  dpi=400,
                  ):
    '''
    311 Plot showing windspeed map and wind speed, along with near sites and transect
    312 relative humidity along transect
    313 water and ice content along transect
    INPUTS:
        cubes need to be passed in
        extentname = { 'waroona' | 'sirivan' } choice of two fire locations
        transect = int from 1 to 6 for transect choice (see figures/transects.png)
        vectorskip reduces quiver density (may need to fiddle)
        quiverscale changes how long the arrows are (also may need to fiddle)
        dtime is datetime of output 
        ext is the plot extension { '.png' | '.eps' }
    '''
    
    # Take required data from cubes (assume already sliced to right time?)
    
    #topog=data['topog']
    #w=data['upward_air_velocity'][tstep]
    #sh=data['specific_humidity'][tstep]
    #s=data['wind_speed'][tstep]
    #u=data['x_wind_destaggered'][tstep]
    #v=data['y_wind_destaggered'][tstep]
    #z=data['zth'][tstep]
    #lat=data['latitude']
    #lon=data['longitude']
    #qc=data['qc'][tstep]
    #Ta=cubes.extract('air_temperature')
    #rh = utils.relative_humidity_from_specific(sh,Ta)
    
    # datetime timestamp for file,title
    dstamp = dtime.strftime("%Y%m%d%H%M")
    stitle = dtime.strftime("%Y %b %d %H:%M (UTC)")
    # figure name and location
    pname="figures/%s/clouds_outline_X%d/fig_%s%s"%(extentname,transect,dstamp,ext)
    
    # set font sizes
    plotting.init_plots()
    
    # get plot extent, and transect
    extent = plotting._extents_[extentname]
    start,end = plotting._transects_["%s%d"%(extentname,transect)]
    
    plt.figure(figsize=[7,10])
    plt.subplot(3,1,1)
    
    # top panel is wind speed surface values
    plotting.map_contourf(extent,s[0],lat,lon,
                          cmap=plotting._cmaps_['windspeed'],
                          clabel='m/s')
    plt.title('Horizontal wind speed')
    

    # start to end x=[lon0,lon1], y=[lat0, lat1]
    plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
             linewidth=2, 
             marker='X', markersize=7,markerfacecolor='white')
    
    # add nearby towns
    if extentname == 'waroona':
        plotting.map_add_locations(['waroona','yarloop'], 
                                   text=['Waroona', 'Yarloop'], 
                                   textcolor='k')
        # add fire ignition
        plotting.map_add_locations(['fire_waroona'],
                                   text = ['Fire ignition'], 
                                   color='r', marker='*', 
                                   textcolor='k')
        # add pyroCB
    else:
        plotting.map_add_locations(['sirivan','uarbry'], 
                                   text=['Sir Ivan','Uarbry'],
                                   dx=[-.02,.05], dy =[-.015,-.03],
                                   textcolor='k')
        # add fire ignition
        plotting.map_add_locations(['fire_sirivan'],
                                   text = ['Fire ignition'], dx=.05,
                                   color='r', marker='*', 
                                   textcolor='k')
        # add pyroCB

    
    # Add vectors for winds
    # just surface, and one every N points to reduce density
    skip = (slice(None,None,vectorskip),slice(None,None,vectorskip))
    #mlon,mlat = np.meshgrid(lon,lat)
    plt.quiver(lon[skip[1]],lat[skip[0]],u[0][skip],v[0][skip], scale=quiverscale)
    
    
    ## Second row is transect plots
    plt.subplot(3,1,2)
    plotting.transect(rh,z,lat,lon,start,end,topog=topog,
                      cmap='plasma',
                      ztop=ztop)
    plt.title("Specific humidity")
    #plotting.transect_w(w,z, lat, lon,start,end,topog=topog)
    
    plt.ylabel('height (m)')
    plt.xlabel('')
    
    plt.subplot(3,1,3)
    plotting.transect_qc(qc,z,lat,lon,start,end,topog=topog,
                        ztop=ztop,)
    # Show transect start and end
    xticks,xlabels = plotting.transect_ticks_labels(start,end)
    plt.xticks(xticks[0::2],xlabels[0::2]) # show transect start and end
    
    # Save figure into animation folder with numeric identifier
    plt.suptitle(stitle)
    print("INFO: Saving figure:",pname)
    plt.savefig(pname,dpi=dpi)
    plt.close()
    return pname

def outline_waroona(pname='figures/site_outline_waroona.png'):
    
    plotting.init_plots()
    extents = plotting._extents_
    latlons = plotting._latlons_
    
    topog,latt,lont = fio.read_topog()
    
    # Google map image tiles view of synoptic map
    fig,ax,proj=plotting.map_google(extents['waroonas'],
                                    zoom=6,
                                    subplotxyn=[2,1,1],
                                    gridlines=[np.arange(-51,-10,2),
                                               np.arange(100,150,4)])
    plt.title("Waroona synoptic")
    ## Add box around zoomed in area
    ax.add_patch(mpatches.Rectangle(xy=latlons['waroona'][::-1], 
                                    width=.4, 
                                    height=.3,
                                    #facecolor=None,
                                    fill=False,
                                    edgecolor='blue',
                                    linewidth=2,
                                    #linestyle='-',
                                    alpha=0.6, 
                                    transform=ccrs.PlateCarree()))
    ## add text?
    
    ## Add scale
    scaleloc=(0.2,0.05)
    plotting.scale_bar(ax,proj,100, location=scaleloc)
    
    
    ## Look at waroona and yarloop
    _,ax2,_ = plotting.map_google(extents['waroona'],zoom=10,fig=fig,subplotxyn=[2,2,3],draw_gridlines=False)
    plt.title("Fire location")
    
    ## Add scale
    plotting.scale_bar(ax2,proj,10, location=scaleloc)
    
    ## Add contour plot showing topography
    plt.subplot(2,2,4)
    plotting.map_topography(extents['waroona'], topog,latt,lont)
    
    plt.savefig(pname)
    print("FIGURE SAVED: ",pname)
    plt.close()
    

def winds_2panel(data,tstep, 
                 extentname='waroona',
                 transect=1, 
                 vectorskip=9,
                 quiverscale=60,
                 ext='.png'
                ):
    '''
    311 Plot showing contour map and wind speed, along with near sites and transect
    312 plot showing vertical motion along transect
    313 plot showing wind speed along transect
    INPUTS:
        dictionary with topog, vert motion, wind speed, wind direction, zth,lat,lon, time
        timestep: data will have a time dimension
        extentname = { 'waroona' | 'sirivan' } choice of two fire locations
        transect = int from 1 to 6 for transect choice
        vectorskip reduces quiver density (may need to fiddle)
        quiverscale changes how long the arrows are (also may need to fiddle)
        ext is the plot extension { '.png' | '.eps' }
    '''
    topog=data['topog']
    ff=None
    if data['firefront'] is not None:
        ff=data['firefront'][tstep]
    w=data['upward_air_velocity'][tstep]
    s=data['wind_speed'][tstep]
    u=data['x_wind_destaggered'][tstep]
    v=data['y_wind_destaggered'][tstep]
    z=data['zth'][tstep]
    lat=data['latitude']
    lon=data['longitude']
    dtime = utils.date_from_gregorian(data['time'])[tstep]
    
    dstamp = dtime.strftime("%Y%m%d%H%M")
    stitle = dtime.strftime("%Y %b %d %H:%M (UTC)")
    
    # figure name and location
    pname="figures/%s/winds_outline_X%d/fig_%s%s"%(extentname,transect,dstamp,ext)
    
    # set font sizes
    plotting.init_plots()
    
    # get plot extent, and transect
    extent = plotting._extents_[extentname]
    start,end = plotting._transects_["%s%d"%(extentname,transect)]
    
    plt.figure(figsize=[7,10])
    plt.subplot(3,1,1)
    
    # top panel is topography
    plotting.map_topography(extent,topog,lat,lon)
    plt.title('Topography, winds')
    
    # Add fire front contour
    if ff is not None:
        plt.contour(lon,lat,np.transpose(ff),np.array([0]), 
                    colors='red')
    
    # start to end x=[lon0,lon1], y=[lat0, lat1]
    plt.plot([start[1],end[1]],[start[0],end[0], ], '--k', 
             linewidth=2, 
             marker='X', markersize=7,markerfacecolor='white')
    
    # add nearby towns
    textcolor='k'
    if extentname == 'waroona':
        plotting.map_add_locations(['waroona','yarloop'], 
                                   text=['Waroona', 'Yarloop'], 
                                   textcolor=textcolor)
        # add fire ignition
        plotting.map_add_locations(['fire_waroona'],
                                   text = ['Fire ignition'], 
                                   color='r', marker='*', 
                                   textcolor=textcolor)
        # add pyroCB
    else:
        plotting.map_add_locations(['sirivan','uarbry'], 
                                   text=['Sir Ivan','Uarbry'],
                                   dx=[-.02,.05], dy =[-.015,-.03],
                                   textcolor=textcolor)
        # add fire ignition
        plotting.map_add_locations(['fire_sirivan'],
                                   text = ['Fire ignition'], dx=.05,
                                   color='r', marker='*', 
                                   textcolor=textcolor)
        # add pyroCB

    
    # Add vectors for winds
    # just surface, and one every 10 points to reduce density
    skip = (slice(None,None,vectorskip),slice(None,None,vectorskip))
    ## colour the arrows
    # map wind speed to colour map domain [0, 1]
    norm = col.Normalize()
    norm.autoscale(s[skip])
    cmap = plt.get_cmap(plotting._cmaps_['windspeed'])
    plt.quiver(lon[skip[1]],lat[skip[0]],u[0][skip],v[0][skip], scale=quiverscale)
               #color=cmap(norm(s[skip])), 
    
    
    ## Second row is transect plots
    plt.subplot(3,1,2)
    plotting.transect_w(w,z, lat, lon,start,end,topog=topog)
    plt.ylabel('height (m)')
    #plt.xlabel('transect')
    
    ax3 = plt.subplot(3,1,3)
    trs, trx, trz = plotting.transect_s(s,z,lat,lon,start,end,topog=topog)
    xticks,xlabels = plotting.transect_ticks_labels(start,end)
    plt.xticks(xticks[0::2],xlabels[0::2]) # show transect start and end
    
    # Annotate max wind speed
    # only care about lower troposphere
    # 60 levels is about 2800m, 80 levels is about 4700m, 70 levels : 3700m
    upto = 70 
    mloc = np.unravel_index(np.argmax(trs[:upto,:],axis=None),trs[:upto,:].shape)
    note="max = %5.1f m/s\n  (at %4.0f m) "%(trs[:upto,:][mloc], trz[:upto,:][mloc])
    trans = ax3.get_xaxis_transform() # x in data untis, y in axes fraction
    ax3.annotate(note, fontsize=15,
                 xy=(0.33, -0.2 ), xycoords=trans)
    
    # Save figure into animation folder with numeric identifier
    plt.suptitle(stitle)
    print("INFO: Saving figure:",pname)
    plt.savefig(pname)
    plt.close()
    return pname
    
    
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
    
    outline_waroona()
    jeff_transect()