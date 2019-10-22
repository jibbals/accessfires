#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 21:04:13 2019

    Try to get fortran PFT calc going
@author: jesse
"""

import numpy as np
from datetime import datetime
import iris
from time import sleep

import matplotlib.pyplot as plt

from utilities import fio, utils, plotting

## from command line I ran python -m numpy.f2py -c HFDiag_calc.f90 -m HFcalc
##     doesn't work on windows (no fortran compilers installed)
## actually I ran f2py3 -c HFDiag_calc.f90 -m HFcalc
## then I renamed the .so file to just HFcalc.so
# this created a .so file callable here
import utilities.fortran.HFcalc as fortranbit
#import HFcalc

## CONSTANTS
phi = 6.0e-5        # Ratio of fire moisture to heat (kg/kg/K)
DbSP = 0.001        # Beta increment along the saturation point curve
ni = 1001           # Number of increments along the SP curve
nj = 20             # Number of iterations for solving plume centre-line height
beta_e = 0.4        # Briggs entrainment paramter (Briggs uses 0.6, we may use up to 0.9)
zp = 1.0            # Number of elements in the pressure (vertical) dimension
Pmin = 35000.0      # Minimum pressure the plume must rise to (Pa)
Tmin = 253.15        # Minimum temperature level the plume must rise to (K).
PorT = 'T'           # Set to 'T' if Tmin is to be used, otherwise set to 'P'
Wmin = 0.05         # Minimum value of mixed-layer wind speed to be considered to influence zb
Umin = 1.0          # Minimum value of mixed-layer horizontal wind speed
betaSPmin = 0.002   # A minimum value of beta to be added to betaSP, to account 
Cpd=1005.7           # Specific heat of dry air (J/kg/K)
Rd=287.04            # Gas constant for dry air (J/kg/K)
LV=2.501e6           # Latent heat of vapourisation (J/kg)
p0=1.0e5             # Standard pressure (Pa)
grav=9.8             # Acceleration due to gravity (m/s^2)
epsln=0.622          # Ratio of gas constants of dry air and water
pi=3.1415926536     # pi
Prcntg = 7          # 1 = 0.0 %, 2 = 5%, 3 = 10% and so on until 21 = 100% (7 used for all experiments up to 23-05-2019)
Dp = 60.0*100.0      # Pressure layer above the surface for which the HDW index is calculated (Pa)
                     # Srock et al. use 500m, ~60 hPa at sea level.  Also grid spacing is 50 hPa above 900 hPa
nlvl = 140             # Number of elements in the pressure (vertical) dimension


'''

! Declare output variables
   real, intent(out) :: UML        ! Mixed-layer horizontal velocity magnitude
   real, intent(out) :: Um         ! Mixed-layer zonal wind
   real, intent(out) :: Vm         ! Mixed-layer meridional wind
   real, intent(out) :: betaFC     ! Free convection beta parameter
   real, intent(out) :: DthFC      ! Free-convection plume excess potential temperature (K)
   real, intent(out) :: zFC        ! Free convection height above ground level (m)
   real, intent(out) :: pFC        ! Free convection pressure (Pa)
   real, intent(out) :: Bflux      ! Net buoyancy flux required for the plume to reach zFC (m^4.s^-3)
   real, intent(out) :: Hflux      ! Net heat flux required for the plume to reach zFC (J.s^-1 = W)

'''

## First read a profile, and surface information required by fortran subroutine
lat,lon = plotting._latlons_['fire_waroona_upwind']
extent = [lon-.01, lon+.01, lat-.01, lat+.01] # just grab real close to latlon
cubes=fio.read_model_run('waroona_old',datetime(2016,1,6,6), extent=extent,
                         add_dewpoint=True, add_winds=True, add_RH=False,
                         add_topog=True, add_z=True, add_theta=True)

## first interpolate everything to latlon
for i in range(len(cubes)):
    cubes[i] = cubes[i].interpolate([('longitude',lon), ('latitude',lat)],
                                    iris.analysis.Linear())
cubedtimes = utils.dates_from_iris(cubes[0])
print("DEBUG:",cubedtimes[0])
print(cubes)

# extract cubes, remove time dim
TTcube = cubes.extract('air_temperature')[0][0]
qq = cubes.extract('specific_humidity')[0][0].data.data # kg/kg
th = cubes.extract('potential_temperature')[0][0].data.data # K
prcube = cubes.extract('air_pressure')[0][0] # Pa?
pr = prcube.data.data # Pa
TT = TTcube.data.data # K
uu,vv,ww = [cube[0].data.data for cube in cubes.extract(['u','v','upward_air_velocity'])] # m/s
dp = cubes.extract('dewpoint_temperature')[0][0].data.data # K

# surface metrics
# surface values in old run are on different time dimension...!?!
zsfc, psfc, Tsfc = cubes.extract(['surface_altitude', 
                                  'surface_air_pressure', 
                                  'surface_temperature'])
zsfc = float(zsfc.data) # m
if len(psfc.shape) > 0:
    psfc = float(psfc[0].data) # Pa
    Tsfc = float(Tsfc[0].data) # K
else:
    psfc = float(psfc.data)
    Tsfc = float(Tsfc.data)

# Find pressure level at which T = Tmin
if PorT == 'T':
    # get first instance where TT is less than Tmin
    Tmin_indices = TT < Tmin
    Pmin = pr[Tmin_indices][0]

if PorT == 'P':
  print('Plume top defined by pressure minimum: ',Pmin/100.0, 'hPa')
elif PorT == 'T':
  print('Plume top defined by Temperature minumum: ',Tmin - 273.15, 'deg C')
else:
  print('******* WARNING: No plume top temperature or pressure has been set.')
  print('Check Pmin, Tmin and PorT')

#print(fortranbit.subroutines.heat_flux_calc.__doc__)
#

skip = slice(None,None,3) # lets make the vert dim less fine for printing simplicity
[TT,qq,uu,vv,ww,th,pr] = [ arr[skip] for arr in [TT,qq,uu,vv,ww,th,pr]]
# Remove values where pressure is lower than 100 hPa
#sub100 = pr>10000
#[TT,qq,uu,vv,ww,th,pr] = [ arr[sub100] for arr in [TT,qq,uu,vv,ww,th,pr]]
nlvl = qq.shape[0] # number of elements in z dimension

"""
plt.plot(th, pr/100., label='th', color='m')
plt.plot(TT, pr/100., label='TT', color='r')
plt.ylim([1000,100])
plt.xlim([200,500])
plt.yscale('log')
plt.legend()
"""

print("DEBUG: RUNNING FORTRAN SUBROUTINE:", )

print("Outside heatflux subroutine:")
print("zsfc,psfc,Tsfc")
print(zsfc,psfc,Tsfc)
print("TT:",TT[:5],"qq:",qq[:5])
print("uu:",uu[:5],"vv:",vv[:5],"ww:",ww[:5])
print("th:",th[:5],"pr:",pr[:5])
print("phi,DbSP,ni,nj,zp,beta_e,Pmin,betaSPmin,Wmin,Umin,Prcntg")
print(phi,DbSP,ni,nj,zp,beta_e,Pmin,betaSPmin,Wmin,Umin,Prcntg)

print("\n\n")
[UML,Um,Vm,betaFC,DthFC,zFC,pFC,Bflux,Hflux]=fortranbit.subroutines.heat_flux_calc(TT,qq,uu,vv,ww,th,pr,nlvl,zsfc,psfc,Tsfc,\
                                                     phi,DbSP,ni,nj,zp,beta_e,Pmin,betaSPmin,Wmin,Umin,Prcntg)

"""
print("\n *** results from heat_flux_calc ***\n")

print("UML, Um, Vm")
print(UML, Um, Vm)
print("betaFC, DthFC, zFC, pFC")
print(betaFC, DthFC, zFC, pFC)
print("Bflux:",Bflux/1e9,"e9")
print(Hflux/1e9, "GWatts")
print(" *** end results from heat_flux_calc ***\n")
"""


#
#Wrapper for ``heat_flux_calc``.
#
#Parameters
#----------

#tt : input rank-1 array('f') with bounds (f2py_tt_d0) # temperature
#qq : input rank-1 array('f') with bounds (f2py_qq_d0) # 
#uu : input rank-1 array('f') with bounds (f2py_uu_d0)
#vv : input rank-1 array('f') with bounds (f2py_vv_d0)
#ww : input rank-1 array('f') with bounds (f2py_ww_d0)
#th : input rank-1 array('f') with bounds (f2py_th_d0)
#pr : input rank-1 array('f') with bounds (f2py_pr_d0)
#nlvl : input int
#zsfc : input float
#psfc : input float
#tsfc : input float
#phi : input float
#dbsp : input float
#ni : input int
#nj : input int
#zp : input float
#beta_e : input float
#pmin : input float
#betaspmin : input float
#wmin : input float
#umin : input float
#prcntg : input int
#
#Returns
#-------
#uml : float
#um : float
#vm : float
#betafc : float
#zfc : float
#pfc : float
#bflux : float
#hflux : float
