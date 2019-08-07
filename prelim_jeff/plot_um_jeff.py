# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 11:59:29 2015

@author: jdk
"""

import matplotlib as mpl
from   matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
from netCDF4 import Dataset
#import sys

import plot_xsec

#mpl.rcParams['font.size'] = 18.0
#mpl.rcParams['axes.formatter.useoffset'] = False

tx = 19

#datadir = 'C:/Users/Jeff/AAA_Work/BNHCRC_CoupledFA/python/data/'
#datadir = '/g/data/en0/mxp548/access-fire/waroona/run1/accessout/'
datadir = '/g/data/en0/mxp548/access-fire/waroona/run3/accessdata/'
datadir2 = '/g/data/en0/mxp548/access-fire/waroona/run2/accessout/'
#ddhh = '0613'
ddhh2 = '0515'
ddhh = '0608'

xsecname='xxx'
zl = 4000

# ends of xsec
lona = 115.7
lata = -32.8
lonb = 116.3
latb = -32.8
lonc = 115.95
latc = -32.6
lond = 115.85
latd = -33.2

dlat = 0.05
nsec = 5

# Read the model data
#ncfile = Dataset(datadir+'umnsaa_pc201601{:s}00.nc'.format(ddhh),'r')
ncfile = Dataset(datadir+'umnsaa_pc201601{:s}.nc'.format(ddhh),'r')
ncfile2 = Dataset(datadir2+'umnsaa_pa201601{:s}.nc'.format(ddhh2),'r')
#zth  = ncfile.variables['height_theta'][0,:,::-1,:] # flip latitudes for interpolate
#zrho = ncfile.variables['height_rho'  ][0,:,::-1,:]
lat  = ncfile.variables['latitude'   ][:]
lon  = ncfile.variables['longitude'  ][:]
lat1 = ncfile.variables['latitude_0' ][:]
lon1 = ncfile.variables['longitude_0'][:]
#z   = ncfile.variables['hybrid_ht'  ][:]
#z1  = ncfile.variables['hybrid_ht_1'][:]
z   = ncfile.variables['model_level_number'  ][:]
z1  = ncfile.variables['model_level_number_0'][:]
#Ta  = ncfile.variables['air_temperature'][0,:,:,:]
Ta  = ncfile.variables['air_temperature'][:,:]
p   = ncfile.variables['air_pressure'][0,0:70,:,:]
#pmsl = ncfile.variables['air_pressure_at_sea_level'][0,0,:,:]
pmsl = ncfile.variables['air_pressure_at_sea_level'][0,:,:]
#u1  = ncfile.variables['x_wind'][0,:,:,:]
#v1  = ncfile.variables['y_wind'][0,:,:,:]
u1  = ncfile.variables['x_wind'][0,0:70,:,:]
v1  = ncfile.variables['y_wind'][0,0:70,:,:]
q   = ncfile.variables['specific_humidity_0'][0,0:70,:,:]
w   = ncfile.variables['upward_air_velocity'][0,0:70,:,:]
#qc  = ncfile.variables['QCL'][0,0:70,:,:] + ncfile.variables['QCF'][0,0:70,:,:]
qc  = ncfile.variables['mass_fraction_of_cloud_liquid_water_in_air'][0,0:70,:,:] + ncfile.variables['mass_fraction_of_cloud_ice_in_air'][0,0:70,:,:]
ncfile.close()

#nz,ny,nx = Ta.shape
nz,ny,nx = p.shape

# Destagger winds
#u = np.tile(np.nan,(nz,ny,nx))
u = np.tile(np.nan,(nz,ny,nx))
u[:,:,1:] = 0.5*(u1[:,:,1:] + u1[:,:,:-1])
v = 0.5*(v1[:,1::,] + v1[:,:-1,:])
s = np.hypot(u,v)
lonu = lon1
latu = lat
lonv = lon
latv = lat1

# Fudge some height data
zth = -(287*300/9.8)*np.log(p/pmsl[np.newaxis,:,:])
zrho = zth

# Fudge some topography data
lont = lon
latt = lat
topog = np.zeros((ny,nx))

# Read the model topography data
#with Dataset(datadir2+'qrparm.orog.mn.nc','r') as ncfile:
with Dataset(datadir2+'umnsaa_pa201601{:s}.nc'.format(ddhh2),'r') as ncfile:
#    topog = ncfile.variables['surface_altitude'][0,0,:,:]
    topog = ncfile.variables['surface_altitude'][:,:]
#    latt = ncfile.variables['rlat' ][:]
#    lont = ncfile.variables['rlon'][:]
    latt = ncfile.variables['latitude' ][:]
    lont = ncfile.variables['longitude'][:]

theta = Ta*(1e5/p)**(287.05/1004.64)

thcon = np.arange(280,320,2)
cmapth = plt.cm.get_cmap('YlOrRd')
   
wcon = 2.0**np.arange(-2,6)
wcon = np.union1d(np.union1d(wcon,-wcon),np.array([0]))
wcon2 = np.array([0])
cmapw = plt.cm.get_cmap('PiYG')
cmapw.set_over('k')
wnorm = col.SymLogNorm(0.25)

scon = np.arange(0,22,2)
cmaps = plt.cm.get_cmap('YlGnBu')

qccon = np.arange(0.0,2.25,0.25)
cmapqc = plt.cm.get_cmap('YlGnBu')

with PdfPages('um_{:s}.pdf'.format(ddhh)) as pdf:
    for isec in range(-(nsec//2),nsec//2+1):
        dl = dlat*isec
        fig = plt.figure(1)
        plt.clf()
        ax = plt.subplot(2,2,1)
        plot_xsec.plot_xsec(lon,lat,theta,lon,lat,zth,lont,latt,topog,
                            (lona,lata+dl),(lonb,latb+dl),
                            thcon,thcon,ztop=zl,cmap=cmapth)
        plt.title('Potential temperature (K)')
        plt.subplot(2,2,2,sharex=ax,sharey=ax)
        plot_xsec.plot_xsec(lon,lat,w,lon,lat,zth,lont,latt,topog,
                            (lona,lata+dl),(lonb,latb+dl),
                            wcon,wcon2,ztop=zl,cmap=cmapw,norm=wnorm,cbformat=tick.ScalarFormatter())
        plt.title('Vertical velocity (m/s)')
        plt.subplot(2,2,3,sharex=ax,sharey=ax)
        plot_xsec.plot_xsec_spd(lonu,latu,u1,lonv,latv,v1,lon,lat,zrho,lont,latt,topog,
                            (lona,lata+dl),(lonb,latb+dl),
                            scon,scon,ztop=zl,cmap=cmaps)
        plt.title('Wind speed (m/s)')
        plt.subplot(2,2,4,sharex=ax,sharey=ax)
        plot_xsec.plot_xsec(lon,lat,qc*1e3,lon,lat,zrho,lont,latt,topog,
                            (lona,lata+dl),(lonb,latb+dl),
                            qccon,np.array([0.0]),ztop=zl,cmap=cmapqc)
        plt.title('Cloud water/ice')
        
        plt.subplots_adjust(right=0.9,hspace=0.4)
        fig.set_size_inches(15,10)
        pdf.savefig(fig)
    
    for isec in range(-(nsec//2),nsec//2+1):
        dl = dlat*isec
        fig = plt.figure(2)
        plt.clf()
        ax = plt.subplot(2,2,1)
        plot_xsec.plot_xsec(lon,lat,theta,lon,lat,zth,lont,latt,topog,
                            (lonc+dl,latc),(lond+dl,latd),
                            thcon,thcon,ztop=zl,cmap=cmapth)
        plt.title('Potential temperature (K)')
        plt.subplot(2,2,2,sharex=ax,sharey=ax)
        plot_xsec.plot_xsec(lon,lat,w,lon,lat,zth,lont,latt,topog,
                            (lonc+dl,latc),(lond+dl,latd),
                            wcon,wcon2,ztop=zl,cmap=cmapw,norm=wnorm,cbformat=tick.ScalarFormatter())
        plt.title('Vertical velocity (m/s)')
        plt.subplot(2,2,3,sharex=ax,sharey=ax)
        plot_xsec.plot_xsec_spd(lonu,latu,u1,lonv,latv,v1,lon,lat,zrho,lont,latt,topog,
                            (lonc+dl,latc),(lond+dl,latd),
                            scon,scon,ztop=zl,cmap=cmaps)
        plt.title('Wind speed (m/s)')
        plt.subplot(2,2,4,sharex=ax,sharey=ax)
        plot_xsec.plot_xsec(lon,lat,qc*1e3,lon,lat,zrho,lont,latt,topog,
                            (lonc+dl,latc),(lond+dl,latd),
                            qccon,np.array([0.0]),ztop=zl,cmap=cmapqc)
        plt.title('Cloud water/ice')
        
        plt.subplots_adjust(right=0.9,hspace=0.4)
        fig.set_size_inches(15,10)
        pdf.savefig(fig)
    
    #for iz in range(10,80,10):
    for iz in range(10,70,10):
        fig = plt.figure(3)
        plt.clf()
        plt.suptitle('level {:02d} (hybrid {:.1f}m)'.format(iz,z[iz]))
        ax2 = plt.subplot(2,2,1)
        plt.contourf(lon,lat,theta[iz,:,:],thcon,cmap=cmapth)
        plt.colorbar()
        plt.plot([lona,lonb],[lata,latb],'k-',linewidth=2)
        plt.plot([lonc,lond],[latc,latd],'k-',linewidth=2)
        plt.title('Potential temperature (K)')
        plt.subplot(2,2,2,sharex=ax2,sharey=ax2)
        plt.contourf(lon,lat,w[iz,:,:],wcon,cmap=cmapw,norm=wnorm)
        plt.colorbar(format=tick.ScalarFormatter())
        plt.contour(lon,lat,w[iz,:,:],[0.0],colors='k')
        plt.plot([lona,lonb],[lata,latb],'k-',linewidth=2)
        plt.plot([lonc,lond],[latc,latd],'k-',linewidth=2)
        plt.title('Vertical velocity (m/s)')
        plt.subplot(2,2,3,sharex=ax2,sharey=ax2)
        plt.contourf(lon,lat,s[iz,:,:],scon,cmap=cmaps)
        plt.colorbar()
        plt.plot([lona,lonb],[lata,latb],'k-',linewidth=2)
        plt.plot([lonc,lond],[latc,latd],'k-',linewidth=2)
        plt.title('Wind speed (m/s)')
        plt.subplot(2,2,4,sharex=ax2,sharey=ax2)
        plt.contourf(lon,lat,qc[iz,:,:]*1e3,qccon,cmap=cmapqc)
        plt.colorbar()
        plt.plot([lona,lonb],[lata,latb],'k-',linewidth=2)
        plt.plot([lonc,lond],[latc,latd],'k-',linewidth=2)
        plt.title('Cloud water/ice')
        fig.set_size_inches(15,10)
        pdf.savefig(fig)
    
    
