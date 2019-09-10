# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 10:16:31 2014
    HISTORY:
        20190910 - jesse: added moist_lapse method
@author: Jeff
"""

import numpy as np
from scipy import integrate
#import matplotlib.pyplot as plt

# Thermodynamic constants, all from Emanuel Convection 
Cpd = 1005.7  # J/kg/K
Cpv = 1870.0  # J/kg/K
Cvd = 719     # J/kg/K
Cvv = 1410    # J/kg/K
Cl  = 4190.0  # J/kg/K # did have 2500???
Rv  = 461.5   # J/kg/K
Rd  = 287.04  # J/kg/K
eps = Rd/Rv
kappa = Rd/Cpd
Lv0 = 2.501e6 # J/kg
g = 9.80665   # m/s**2

# Cpm ####################################################################
def Cpm(r):
    # Heat capacity at constant pressure of moist air with mixing ratio r
    return (Cpd + Cpv*r)/(1 + r)
    
# Cvm ####################################################################
def Cvm(r):
    # Heat capacity at constant volume of moist air with mixing ratio r
    return (Cvd + Cvv*r)/(1 + r)
    
# Lv #####################################################################
def Lv(T):
    # Latent heat of vapourisation of water
    return Lv0 + (Cpv - Cl)*(T - 273.15)
    
# e_sat ##################################################################
def e_sat(T):
    # Calculates the saturation vapour pressure in Pa given T in K
    TC = T - 273.15
    es = 611.2*np.exp( 17.67*TC/(243.5 + TC) )
    return es

# e_satd ##################################################################
def e_satd(T):
    # Calculates the saturation vapour pressure in Pa given T in K and
    # its derivative with respect to T
    TC = T - 273.15
    es = 611.2*np.exp( 17.67*TC/(243.5 + TC) )
    desdT = 243.5*es / (243.5 + TC)**2
    return es,desdT

# q_sat ##################################################################
def q_sat(T, p=1e5):
# q_sat(T,p) gives temperature T (K) and pressure p (Pa) returns the 
# saturation specific humidity in kg/kg over pure water. 
# q_sat(T) takes the pressure as 1000 hPa.
# Ref: Buck, JAM 20, 1981, 1527-1532.

    esat = e_sat(T)
    qs = eps * esat / (p - (1-eps)*esat)
    return qs

# r_sat ##################################################################
def r_sat(T, p=1e5):
# r_sat(T,p) gives temperature T (K) and pressure p (Pa) returns the 
# saturation mixing ratio in kg/kg over pure water. 
# r_sat(T) takes the pressure as 1000 hPa.
# Ref: Buck, JAM 20, 1981, 1527-1532.

    esat = e_sat(T)
    rsat = eps * esat / (p - esat)
    return rsat

# r_satd ##################################################################
def r_satd(T, p=1e5):
# r_sat(T,p) gives temperature T (K) and pressure p (Pa) returns the 
# saturation mixing ratio in kg/kg over pure water and tis derivative
# with respect to T. 
# r_sat(T) takes the pressure as 1000 hPa.
# Ref: Buck, JAM 20, 1981, 1527-1532.

    esat,desdT = e_satd(T)
    rsat = eps * esat / (p - esat)
    drsdT = desdT*eps*p/(p - esat)**2
    return rsat,drsdT

# theta_e ################################################################
def theta_e(T, r = 0.0, p = 1e5):
    # th_e = theta_e(T, r, p)
    # Returns the equivalent potential temperature, given the temperature (in K),
    # mixing ratio (in kg/kg)  and pressure (in Pa). Uses the formulae of
    # Bolton, MWR, 1980. 

    if np.all(r==0):
        th_e = T * ((1e5/p) ** 0.2854)  
    else:
        es = e_sat(T)
        rs = eps*es/(p - es)
        r = np.minimum(r, rs)
        e = p * r / (eps + r)
        loge = np.log( 0.01*np.maximum(e, 1e-100) )
        T_lcl = 2840.0 / (3.5*np.log(T) - loge - 4.805) + 55.0
        th_e = ( T * ((1e5/p) ** (0.2854 - 0.079912*r)) 
                   * np.exp((3376/T_lcl - 2.54) * r * (1.0 + 0.81*r)) )
        #print(r,e,T_lcl,th_e)
    
    return th_e

# T_d ########################################################################
def T_d(e=None, r=None, p=None):
    # Need either e, or r and p. e and p in Pa
    
    if (e == None):
        if (r == None or p == None):
            print('Error in Td')
            return np.nan
        else:
            e = p*r/(eps + r)
    loge = np.log(e/611.2)
    return 243.5*loge / (17.67 - loge) + 273.15

# T_v ########################################################################
def T_v(T, r):
    return T*(1 + r/eps)/(1 + r)
    
# lift_parcel ################################################################
def lift_parcel(Tp,rp,pp, p):
     # lift the parcel (Tp,rp,pp) to pressures p (in Pa). Return T,r

    n = p.size
    
    # Default values   
    iflag = 0

    # Thermodynamic constant
    CpvMCl = Cpv - Cl

    # Define various parcel quantities, including parcel reversible entropy, Sp                         ***
    Esp = e_sat(Tp)
    Evp = rp*pp/(eps + rp)
    RH = min(Evp/Esp, 1.0)
    Lv = Lv0 - CpvMCl*(Tp - 273.15)
    Sp = (Cpd + rp*Cl)*np.log(Tp) - Rd*np.log(1e-2*(pp-Evp)) + Lv*rp/Tp - rp*Rv*np.log(RH)
    
    # Find lifted condensation pressure, pLCL   ***
    CHI = Tp / (1669.0 - 122.0*RH - Tp)
    pLCL = pp * (RH**CHI)
    Tlcl = 2840.0 / (3.5*np.log(Tp) - np.log(Evp*1e-2) - 4.805) + 55.0 # from Bolton

    # Updraft loop
    Tup = np.zeros(p.shape)
    rup = np.zeros(p.shape)
    plast = pLCL; Tlast = Tlcl  # used to guess parcel temps in updraft loop
    
    for j in np.arange(n):

        if p[j] >= pLCL:
            # Parcel quantities below lifted condensation level
            Tup[j] = Tp*(p[j]/pp)**(Rd/Cpd)
            rup[j] = rp
        else:
            # Parcel quantities above lifted condensation level
            # First, use SALR to get guess for new parcel temperature
            ES = e_sat(Tlast)
            rg = eps*ES/(p[j] - ES)
            RCp = (Rd*(1 + rg/eps))/(Cpd + Cpv*rg)
            Lterm = (Lv*rg)/(Rd*Tlast)
            gam = (Tlast/plast)*RCp*(1 + Lterm)/(1 + Lterm*RCp*Lv/(Rv*Tlast))
            Tg = Tlast - (plast - p[j])*gam
            ES = e_sat(Tg)
            rg = eps*ES/(p[j] - ES)

            # Iteratively calculate lifted parcel temperature and mixing   ***
            # ratio for reversible ascent
            iflag = 2
            for nc in range(500):
                # Calculate estimates of the rates of change of the entroy
                # with temperature at constant pressure
                Lv = Lv0 - CpvMCl*(Tg - 273.15)
                SL = (Cpd + rp*Cl + Lv*Lv*rg/(Rv*Tg*Tg))/Tg # approx dS/dTg
                EM = rg*p[j] / (eps + rg)
                SG = (Cpd + rp*Cl)*np.log(Tg) - Rd*np.log(1e-2*(p[j] - EM)) + Lv*rg/Tg
                TgOLD = Tg
                Tg = Tg + (Sp - SG)/SL  # Newton-Raphson correction
                esg = e_sat(Tg)
                rg = eps*esg / (p[j] - esg)           

                # Test for convergence
                if abs(TgOLD - Tg) < 1e-4:
                    iflag = 0
                    break
            Tlast = Tg
            plast = p[j]
                                
            Tup[j] = Tg
            rup[j] = rg

            if iflag == 2:
                print('Convergence failure in lift_parcel',nc,Tg)
                return Tup,rup,iflag
 
    return Tup,rup,iflag

# CAPE #######################################################################
def CAPE(Tp,rp,pp, T,r,p):
    # cape,Tnb,iflag = cape(Tp,rp,pp, T,r,p)

    # Calculates the cape of a parcel.
    # Inputs:
    # Parcel properties are pressure pp (mb), temperature Tp (K) and 
    # mixing ratio rp (kg/kg).
    # Environmental sounding temperature T (K) and mixing ratio r (kg/kg) 
    # and pressure (p in mb). 
    # Outputs:
    # cape is cape 
    # Tnb is the temperature at the level of neutral buoyancy
    # iflag is a flag:
    #   iflag = 0, success
    #   iflag = 1, bad input data
    #   iflag = 2, parcel temperature did not converge in moist ascent
    # Based on Kerry Emanuel's fortran code

    n = T.size
    
    # Default values   
    cape = 0.0
    Tnb = T[0]
    iflag = 1

    # Thermodynamic constant
    CpvMCl = Cpv - Cl

    # Sanity check
    if (rp < 1e-6 or Tp < 200.0 or pp > p[0]):
       iflag = 1
       return cape,Tnb,iflag

    # Define various parcel quantities, including parcel reversible entropy, Sp                         ***
    Esp = 1e-2*e_sat(Tp)
    Evp = rp*pp/(eps + rp)
    RH = min(Evp/Esp, 1.0)
    Lv = Lv0 - CpvMCl*(Tp - 273.15)
    Sp = (Cpd + rp*Cl)*np.log(Tp) - Rd*np.log(pp-Evp) + Lv*rp/Tp - rp*Rv*np.log(RH)
    
    # Find lifted condensation pressure, pLCL   ***
    CHI = Tp / (1669.0 - 122.0*RH - Tp)
    pLCL = pp * (RH**CHI)
    Tlcl = 2840.0 / (3.5*np.log(Tp) - np.log(Evp) - 4.805) + 55.0 # from Bolton

    # Updraft loop
    delTv = np.zeros(p.shape)
    jMIN = 10000
    Tup = np.zeros(p.shape)
    rup = np.zeros(p.shape)
    plast = pLCL; Tlast = Tlcl  # used to guess parcel temps in updraft loop
    
    for j in np.arange(n):

        # Don't bother lifting parcel above 60 mb and skip sections of sounding below parcel level  ***
        if p[j] < 59.0 or p[j] > pp:
            continue

        jMIN = min(jMIN,j)

        if p[j] >= pLCL:
            # Parcel quantities below lifted condensation level
            Tg = Tp*(p[j]/pp)**(Rd/Cpd)
            rg = rp
            Tup[j] = Tg
            rup[j] = rg

            # Calculate buoyancy
            TLVR = Tg*(1 + rg/eps)/(1 + rg)
            delTv[j] = TLVR - T[j]*(1 + r[j]/eps)/(1 + r[j])
        else:
            # Parcel quantities above lifted condensation level
            # First, use SALR to get guess for new parcel temperature
            ES = 1e-2*e_sat(T[j]) # 1e-2 to convert to hPa, same as p
            rg = eps*ES/(p[j] - ES)
            RCp = (Rd*(1 + rg/eps))/(Cpd + Cpv*rg)
            Lterm = (Lv*rg)/(Rd*Tlast)
            gam = (Tlast/plast)*RCp*(1 + Lterm)/(1 + Lterm*RCp*Lv/(Rv*Tlast))
            Tg = Tlast - (plast - p[j])*gam
            #Tg = T[j]          
            ES = 1e-2*e_sat(Tg) # 1e-2 to convert to hPa, same as p
            rg = eps*ES/(p[j] - ES)

            # Iteratively calculate lifted parcel temperature and mixing   ***
            # ratio for reversible ascent
            iflag = 2
            for nc in range(500):
                # Calculate estimates of the rates of change of the entroy
                # with temperature at constant pressure
                Lv = Lv0 - CpvMCl*(Tg - 273.15)
                SL = (Cpd + rp*Cl + Lv*Lv*rg/(Rv*Tg*Tg))/Tg # approx dS/dTg
                EM = rg*p[j] / (eps + rg)
                SG = (Cpd + rp*Cl)*np.log(Tg) - Rd*np.log(p[j] - EM) + Lv*rg/Tg
                #if nc < 3: Remove step limiter in the hope that a better guess will fix
                #    AP=0.3
                #else:
                #    AP=1.0
                TgOLD = Tg
                #Tg = Tg + AP*(Sp - SG)/SL  # Newton-Raphson correction
                Tg = Tg + (Sp - SG)/SL  # Newton-Raphson correction
                esg = 1e-2*e_sat(Tg)
                rg = eps*esg / (p[j] - esg)           

                # Test for convergence
                if abs(TgOLD - Tg) < 1e-4:
                    iflag = 0
                    break
            Tlast = Tg
            plast = p[j]
                                
            Tup[j] = Tg
            rup[j] = rg

            if iflag == 2:
                print('Convergence failure in cape routine',nc,Tg)
                return cape,Tnb,iflag
 
            # Calculate buoyancy - uncomment one of these two lines
            TLVR = Tg*(1 + rg/eps)/(1 + rp)   # reversible
            # TLVR = Tg*(1 + rg/eps)/(1 + rg)   # pseudo-adiabatic
            delTv[j] = TLVR - T[j]*(1 + r[j]/eps) / (1 + r[j])

    # Begin loop to find NA, PA, and cape from reversible ascent ***
    NA = 0.0
    PA = 0.0

    # Find maximum level of positive buoyancy, INB
    INB = 0
    for j in range(jMIN,n):
        if delTv[j] > 0.0:
            INB = max(INB,j)

    #print('jmin>>',jMIN,INB)
    # Find positive and negative areas and cape
    if INB > 0:
        # Find area between parcel pressure and first level above it
        # Changed from Emanuel's version to  allow the parcel to be different
        # from the environment, Emanuel assumed they were the same.
        ## PFAC = Rd * (pp - p[jMIN]) / (0.5*(pp + p[jMIN]))
        ## PA = PA + PFAC*max( delTv[jMIN], 0.0 )
        ## NA = NA - PFAC*min( delTv[jMIN], 0.0 )
        # Environment interpolated to parcel pressure
        if pp > p[jMIN]:
            if jMIN == 0:
                T_pp = T[jMIN]
                r_pp = r[jMIN]
            else:
                wt = (pp - p[jMIN])/(p[jMIN-1] - p[jMIN])
                T_pp = wt*T[jMIN-1] + (1-wt)*T[jMIN]
                r_pp = wt*r[jMIN-1] + (1-wt)*r[jMIN]
            Tv_pp = T_pp*(1 + r_pp/eps)/(1 + r_pp) # environment Tv at pp
            Tvp   = Tp  *(1 + rp  /eps)/(1 + rp  ) # parcel Tv
            PFAC = Rd*(Tvp - Tv_pp + delTv[jMIN]) * (pp - p[jMIN]) / (pp + p[jMIN])
            PA = PA + max(PFAC,0.0)
            NA = NA - min(PFAC,0.0)
        #print('initial bit',PA,NA)

        for j in range(jMIN+1,INB+1):
            PFAC = Rd*(delTv[j] + delTv[j-1]) * (p[j-1] - p[j]) / (p[j] + p[j-1])
            PA = PA + max(PFAC,0.0)
            NA = NA - min(PFAC,0.0)
            #print('int loop',j,PA,NA)
  
        # Find residual positive area above INB and TO
        PAT = 0.0
        Tnb = T[INB]
        if INB < n-1:
            PINB = (p[INB+1]*delTv[INB] - p[INB]*delTv[INB+1])/(delTv[INB] - delTv[INB+1])
            PAT = Rd*delTv[INB]*(p[INB] - PINB)/(p[INB] + PINB)
            Tnb = (T[INB]*(PINB - p[INB+1]) + T[INB+1]*(p[INB] - PINB))/(p[INB] - p[INB+1])

        # Find cape
        #print('final bit',PA,PAT,NA)
        cape = max( PA+PAT-NA, 0.0 )

    return cape, Tnb, iflag#, Tup
    
"""
Created on Tue Sep 10 13:22:58 2019
    moist_lapse rate modified from https://unidata.github.io/MetPy/latest/_modules/metpy/calc/thermo.html
    this seems to be something like moist equivalent potential temp
@author: jgreensl
"""
def moist_lapse(pressure, temperature, ref_pressure=None):
    r"""Calculate the temperature at a level assuming liquid saturation processes.
    
    Pressures in Pa, temp in Kelvin
    
    This function lifts a parcel starting at `temperature`. The starting pressure can
    be given by `ref_pressure`. Essentially, this function is calculating moist
    pseudo-adiabats.

    Parameters
    ----------
    pressure : `pint.Quantity`
        The atmospheric pressure level(s) of interest
    temperature : `pint.Quantity`
        The starting temperature
    ref_pressure : `pint.Quantity`, optional
        The reference pressure. If not given, it defaults to the first element of the
        pressure array.

    Returns
    -------
    `pint.Quantity`
       The temperature corresponding to the starting temperature and
       pressure levels.

    See Also
    --------
    dry_lapse : Calculate parcel temperature assuming dry adiabatic processes
    parcel_profile : Calculate complete parcel profile

    Notes
    -----
    This function is implemented by integrating the following differential
    equation:

    .. math:: \frac{dT}{dP} = \frac{1}{P} \frac{R_d T + L_v r_s}
                                {C_{pd} + \frac{L_v^2 r_s \epsilon}{R_d T^2}}

    This equation comes from [Bakhshaii2013]_.

    """
    def dt(t, p):
        # make sure t is in (K) and p is in (Pa)
        #t = units.Quantity(t, temperature.units)
        #p = units.Quantity(p, pressure.units)
        rs = r_sat(t,p) # T (K) and pressure p (Pa) returns the saturation mixing ratio in kg/kg
        #frac = ((mpconsts.Rd * t + mpconsts.Lv * rs)
        #        / (mpconsts.Cp_d + (mpconsts.Lv * mpconsts.Lv * rs * mpconsts.epsilon
        #                            / (mpconsts.Rd * t * t)))).to('kelvin')
        frac = ((Rd * t + Lv0 * rs)
                / (Cpd + (Lv0 * Lv0 * rs * eps
                                    / (Rd * t * t))))
        #JESSE TODO: Make sure frac is in Kelvin!!
        return frac / p

    if ref_pressure is None:
        ref_pressure = pressure[0]

    #pressure = pressure.to('mbar')
    pressure = pressure/100. # Pa to hPa
    #ref_pressure = ref_pressure.to('mbar')
    ref_pressure = ref_pressure/100. # Pa to hPa
    #temperature = atleast_1d(temperature)

    side = 'left'

    pres_decreasing = (pressure[0] > pressure[-1])
    if pres_decreasing:
        # Everything is easier if pressures are in increasing order
        pressure = pressure[::-1]
        side = 'right'

    ref_pres_idx = np.searchsorted(pressure.m, ref_pressure.m, side=side)

    ret_temperatures = np.empty((0, temperature.shape[0]))

    if ref_pressure > pressure.min():
        # Integrate downward in pressure
        pres_down = np.append(ref_pressure, pressure[(ref_pres_idx - 1)::-1])
        trace_down = integrate.odeint(dt, temperature.squeeze(), pres_down.squeeze())
        ret_temperatures = np.concatenate((ret_temperatures, trace_down[:0:-1]))

    if ref_pressure < pressure.max():
        # Integrate upward in pressure
        pres_up = np.append(ref_pressure, pressure[ref_pres_idx:])
        trace_up = integrate.odeint(dt, temperature.squeeze(), pres_up.squeeze())
        ret_temperatures = np.concatenate((ret_temperatures, trace_up[1:]))

    if pres_decreasing:
        ret_temperatures = ret_temperatures[::-1]

    return ret_temperatures.T.squeeze()

