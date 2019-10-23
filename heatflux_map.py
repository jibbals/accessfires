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
from platform import platform # check for windows
from metpy.plots import SkewT


import matplotlib.pyplot as plt

from utilities import fio, utils, plotting

## from command line I ran python -m numpy.f2py -c HFDiag_calc.f90 -m HFcalc
##     doesn't work on windows (no fortran compilers installed)
## actually I ran f2py3 -c HFDiag_calc.f90 -m HFcalc
## then I renamed the .so file to just HFcalc.so
# this created a .so file callable here

onlaptop=True
if 'Windows' in platform():
    onlaptop=False
else:
    import utilities.fortran.HFcalc as fortranbit
import utilities.fortran.PFT as HFj

    
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

# extract cubes, remove time dim
TTcube = cubes.extract('air_temperature')[0][0]
qq = cubes.extract('specific_humidity')[0][0].data.data # kg/kg
th = cubes.extract('potential_temperature')[0][0].data.data # K
prcube = cubes.extract('air_pressure')[0][0] # Pa
pr = prcube.data.data # Pa
TT = TTcube.data.data # K
Td = cubes.extract('dewpoint_temperature')[0][0].data.data # K
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
[TT,qq,uu,vv,ww,th,pr,Td] = [ arr[skip] for arr in [TT,qq,uu,vv,ww,th,pr,Td]]
# Remove values where pressure is lower than 100 hPa
#sub100 = pr>10000
#[TT,qq,uu,vv,ww,th,pr] = [ arr[sub100] for arr in [TT,qq,uu,vv,ww,th,pr]]
nlvl = qq.shape[0] # number of elements in z dimension


[UMLp,Ump,Vmp,betaFCp,DthFCp,zFCp,pFCp,Bfluxp,Hfluxp]=HFj.PFT(TT,qq,uu,vv,ww,th,pr, zsfc,psfc,Tsfc, phi,DbSP,ni,nj,zp,beta_e,Pmin,betaSPmin,Wmin,Umin,Prcntg-1)

#### Run also on fortran code to compare..
if onlaptop:
    
    
    print("DEBUG: RUNNING FORTRAN SUBROUTINE:", )

    print("\n\n")
    [UML,Um,Vm,betaFC,DthFC,zFC,pFC,Bflux,Hflux]=fortranbit.subroutines.heat_flux_calc(TT,qq,uu,vv,ww,th,pr,nlvl,zsfc,psfc,Tsfc,\
                                                     phi,DbSP,ni,nj,zp,beta_e,Pmin,betaSPmin,Wmin,Umin,Prcntg)
    
    """
    print("\n *** results from fortran heat_flux_calc ***\n")
    
    print("UML, Um, Vm")
    print(UML, Um, Vm)
    print("betaFC, DthFC, zFC, pFC")
    print(betaFC, DthFC, zFC, pFC)
    print("Bflux:",Bflux/1e9,"e9")
    print(Hflux/1e9, "GWatts")
    print(" *** end results from heat_flux_calc ***\n")

print("\n *** results from PORTED heat_flux_calc ***\n")    
print("UML, Um, Vm")
print(UMLp, Ump, Vmp)
print("betaFC, DthFC, zFC, pFC")
print(betaFCp, DthFCp, zFCp, pFCp)
print("Bflux:",Bfluxp/1e9,"e9")
print(Hfluxp/1e9, "GWatts")
print(" *** end results from PORTED heat_flux_calc ***\n")
"""

skew = SkewT(rotation=45)
skew.plot(pr/100.,TT-273.15,'r', linewidth=2)
skew.plot(pr/100.,Td-273.15,'m', linewidth=2)
#skew.plot(pr/100.,th-273.15,'k', linewidth=2)
skew.plot_moist_adiabats()
#skew.plot_dry_adiabats()
# set limits
skew.ax.set_ylim(1000,100)
skew.ax.set_xlim(-30,60)
plt.show()