# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 13:32:37 2019
    Look at wind vector fields in 3D
@author: jgreensl
"""

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.ticker as tick
import numpy as np
from datetime import datetime

from utilities import fio, plotting, constants

def transect_vector_projection():
    '''
    Trying to project wind vectors onto transect plane
    TODO: Not yet working I think
    '''
    # get plot extent, and transect
    extentname='waroona'
    extent = constants.extents[extentname]
    transect=3
    start,end = plotting._transects_["%s%d"%(extentname,transect)]
    
    ## Read topography, and winds
    tstep = 0
    slv,ro1,th1,th2 = fio.read_waroona(dtime0, extent=extent, add_winds=True)
    z, = ro1.extract('z_ro')
    topog, = slv.extract('topog')
    u,v,s = ro1.extract(['u','v','s'])
    qc, = th2.extract('qc')
    w, = th1.extract('upward_air_velocity')
    
    
    # Get sign and magnitude of u,v along x axis (transect vector)
    # u is along lons, v is along lats, start,end are [lat,lon]
    ui = u[tstep].data
    vi = v[tstep].data
    wi = w[tstep].data
    lat,lon = u.coord('latitude').points,u.coord('longitude').points
    zi = z.data
    n_points_X = 40
    
    ut = utils.cross_section(ui,lat,lon,start,end,npoints=n_points_X)
    vt = utils.cross_section(vi,lat,lon,start,end,npoints=n_points_X)
    Xprojt = (vi*(end[1]-start[1]) + ui*(end[0]-start[0])) / np.sqrt((end[0]-start[0])**2+(end[1]-start[1])**2)
    Xproj0 = Xprojt * (end[0]-start[0])
    Xproj1 = Xprojt * (end[1]-start[1])
    X=np.sqrt(Xproj0**2 + Xproj1**2)  # horizontal magnitude along slice dim
    negs=Xproj0+Xproj1 < 0.0 # where is it negative
    X[negs] = -X[negs]
    tranX = utils.cross_section(X,lat,lon,start=start,end=end,npoints=n_points_X)
    tranZ = utils.cross_section(wi,lat,lon,start=start,end=end,npoints=n_points_X)
    
    # show w transect
    plt.close()
    
    # streamplot doesn't work with non-grid inputs (Y changes as X changes...)
    # can plot using model level height approximation
    #mh = np.nanmean(slz,axis=1)
    #slmh = np.zeros(slz.shape) + mh[:,np.newaxis]
    #plt.streamplot(slx[skip2d], slmh[skip2d], tranX[skip2d], tranZ[skip2d])
    
    plt.figure(figsize=[7,10])
    plt.subplot(3,1,1)
    
    # top panel is topography
    plotting.map_topography(extent,topog,lat,lon)
    plt.title('Topography, winds')
    
    # Add fire front contour
    with warnings.catch_warnings():
        # ignore warning when there are no fires:
        warnings.simplefilter('ignore')
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
    plt.get_cmap(plotting._cmaps_['windspeed'])
    plt.quiver(lon[skip[1]],lat[skip[0]],u[0][skip],v[0][skip], scale=quiverscale)
               #color=cmap(norm(s[skip])), 
    
    
    ## Second row is transect plots
    plt.subplot(3,1,2)
    wslice, xslice, zslice  = plotting.transect_w(w,z, lat, lon,start,end,
                                                  npoints=100,topog=topog,
                                                  lines=None)
    plt.ylabel('height (m)')
    ## Add contour where clouds occur
    qcslice = utils.cross_section(qc,lat,lon,start,end, npoints=100)
    with warnings.catch_warnings():
        # ignore warning when there are no clouds:
        warnings.simplefilter('ignore')
    plt.contour(xslice,zslice,qcslice,np.array([cloud_threshold]),colors='teal')
    
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
    
    
# Read cloud and wind data
def try_3d():
    '''
    Show topography and wind vectors in a 3d framework
    '''
    slv, ro1, th1, th2 = fio.read_waroona(dtime= datetime(2016,1,6,9),extent=constants.extents['waroona'], add_winds=True)
    
    topog, = slv.extract(['topog'])
    u,v = ro1.extract(['u','v'])
    w   = th1.extract(['upward_air_velocity'])[0]
    zth = th1.extract('z_th')[0] 
    # show 1 in 10 lats and lons, and just show specific heights
    skip = (slice(None,None,10), slice(None,None,10), np.array([10,28,40,50]))
    
    lon = u.coord('longitude').points
    lat = u.coord('latitude').points
    lh = u.coord('level_height').points # approximate level height
    mlon2,mlat2 = np.meshgrid(lon,lat)
    mlon3,mlat3,mlh3 = np.meshgrid(lon,lat,lh)
    #for thing in [u,v,w,zth,lon,lat,lh,mlon,mlat,mlh]:
    #    print(np.shape(thing))
    
    # Need to rotate so z dim is last:
    u3 = np.moveaxis(u[0].data, 0,2)
    v3 = np.moveaxis(v[0].data, 0,2)
    w3 = np.moveaxis(w[0].data, 0,2)
    zth3=np.moveaxis(zth.data,0,2)
    mlon3s, mlat3s, mlh3s = mlon3[skip],mlat3[skip],mlh3[skip]
    u3s = u3[skip]
    v3s = v3[skip]
    w3s = w3[skip]
    zth3s = zth3[skip]
    
    # w3 will form the colour map
    # flatten and normalize
    #w3n = (w3s.ravel() - w3s.min()) / w3s.ptp()
    w3sr = w3s.ravel()
    # repeat for each body line and two head lines
    w3src = np.concatenate((w3sr, np.repeat(w3sr,2)))
    
    plt.close() 
    fig=plt.figure(figsize=[10,10])
    #ax=plt.subplot(1,2,1, projection='3d')
    ax = fig.gca(projection='3d')
    
    # Plot topography surface
    ax.plot_surface(mlon2, mlat2, topog.data, cmap='terrain', 
                    linewidth=0, antialiased=True, vmin=-150,vmax=550)
    
    # Show wind vectors overhead, coloured by vertical velocity
    quiv= plt.quiver(mlon3s,mlat3s,zth3s, u3s, v3s, w3s,
                     cmap=plotting._cmaps_['verticalvelocity'],
                     norm=col.SymLogNorm(0.5),
                     length=0.025,normalize=True)
    # this seems to set the colours!
    quiv.set_array(w3src)
    cb=plt.colorbar(format=tick.ScalarFormatter())
    cb.set_label('vertical motion (m/s)')
    
    # lower camera elevation
    ax.view_init(elev=14)
    plt.savefig('figures/waroona/wind_vectors1.png')
    #plt.close()




