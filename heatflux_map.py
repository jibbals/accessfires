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

from utilities import fio, utils, plotting

## from command line I ran python -m numpy.f2py -c HFDiag_calc.f90 -m HFcalc
##     doesn't work on windows (no fortran compilers installed)
# this created a .so file callable here
import utilities.fortran.HFcalc as fortranbit

'''
## Declare additional single-column arrays
# lvl is vertical levels array
qq = lvl*0.0               # Specific Humidity. Units kg/kg
th = lvl*0.0               # Potential temperature. Units K
pr = lvl*100.0             # Pressure. Units Pa
TT = lvl*0.0               # Air Temperature. Units K
uu = lvl*0.0               # Zonal wind. Units m/s
vv = lvl*0.0               # Meridional wind. Units m/s
ww = lvl*0.0               # Vertical wind. Units m/s
dp = lvl*0.0               # Dew point temperature.  Units K
'''

## First read a profile, and surface information required by fortran subroutine
lat,lon = plotting._latlons_['upwind']
extent = [lon-.02, lon+.02, lat-.02, lat+.02] # just grab real close to latlon
cubes=fio.read_model_run('waroona_old',datetime(2016,1,6,5), extent=extent,
                         add_dewpoint=True, add_winds=True, add_RH=False,
                         add_topog=True, add_z=True)

## first interpolate everything to latlon
for i in range(len(cubes)):
    cubes[i] = cubes[i].interpolate([('longitude',lon), ('latitude',lat)],
                                    iris.analysis.Linear())
print("debug:",cubes)


#press.convert_units('hPa') # convert to hPa
#Temp.convert_units('K') # convert to sir Kelvin
#Tempd.convert_units('K') # convert to kelvin
    
#
## Plot T, and Td
## pull out data array (units don't work with masked arrays)
#p = np.squeeze(press0.data.data) * units(str(press.units))
#p = p.to(units.mbar)
#T = np.squeeze(temp0.data.data) * units(str(Temp.units))
#T = T.to(units.degC)
#Td = np.squeeze(tempd0.data.data) * units(str(Tempd.units))
#Td = Td.to(units.degC)
    
#heat_flux_calc(TT,qq,uu,vv,ww,th,pr,pd,zsfc,psfc,Tsfc, &
#               phi,DbSP,ni,nj,zp,beta_e,Pmin,betaSPmin,Wmin,Umin,Prcntg, &
#               UML,Um,Vm,betaFC,zFC,pFC,Bflux,Hflux)



#
print(fortranbit.subroutines.heat_flux_calc.__doc__)
#
#uml,um,vm,betafc,zfc,pfc,bflux,hflux = heat_flux_calc(tt,qq,uu,vv,ww,th,pr,pd,zsfc,psfc,tsfc,phi,dbsp,ni,nj,zp,beta_e,pmin,betaspmin,wmin,umin,prcntg)
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
#pd : input int
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