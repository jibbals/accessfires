#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 10:00:09 2019

    Copy kevin's code but into python

@author: jesse
"""

###
## IMPORTS
###
import numpy


###
## GLOBALS
###

Cpd=1005.7           # Specific heat of dry air (J/kg/K)
Rd=287.04            # Gas constant for dry air (J/kg/K)
kappa=Rd/Cpd         # Ratio of gas constant to specific heat
LV=2.501e6           # Latent heat of vapourisation (J/kg)
P0=1.0e5             # Standard pressure (Pa)
grav=9.8             # Acceleration due to gravity (m/s^2)
epsln=0.622          # Ratio of gas constants of dry air and water
pi=3.1415926536      # pi


###
## METHODS
###

def PFT(TT,qq,uu,vv,ww,th,pr,
        zsfc,psfc,Tsfc,
        phi,DbSP,ni,nj,zp,beta_e,Pmin,betaSPmin,Wmin,Umin,
        Prcntg=6):
    """Calculate required energy to form a pyroCB.
    
    Port of Kevin Tory fortran code to calculate PFT
    
    TODO: tidy up these comments
       integer, intent(in) :: pd       ! Number of elements in the pressure (vertical) dimension
       real, intent(in) :: zsfc        ! Surface height (m)
       real, intent(in) :: psfc        ! Surface pressure (Pa)
       real, intent(in) :: Tsfc        ! Surface temperature (K)
       real, intent(in) :: phi         ! Ratio of fire moisture to heat (kg/kg/K)
       real, intent(in) :: DbSP        ! Beta increment along the saturation point curve
       integer, intent(in) :: ni       ! Number of increments along the SP curve
       integer, intent(in) :: nj       ! Number of iterations for solving plume centre-line height
       real, intent(in) :: beta_e      ! Briggs entrainment paramter (Briggs uses 0.6, we may use up to 0.9)
       real, intent(in) :: zp          ! The power to which zz is raised in calculating wt. (2/3 for MTT, 1.0 for Briggs)
       real, intent(in) :: Pmin        ! The minimum pressure the moist plume needs to rise to for 
                                       !   pyroCb to be considered.  Probably choose 300 to 250 hPa, 
                                       !   just below the tropopause.
       real, intent(in) :: betaSPmin   ! A minimum value of beta to be added to betaSP, to account 
                                       !   for buoyancy losses from entrainment etc. above the condensation level
       real, intent(in) :: Wmin        ! Minimum value of mixed-layer wind speed to be considered to influence zb
       real, intent(in) :: Umin        ! Minimum value of mixed-layer horizontal wind speed
       integer, intent(in) :: Prcntg   ! Index for specifying percentage of plume required to reach zFC
     
       real, intent(in), dimension(:) :: TT       ! Column temperature (K)
       real, intent(in), dimension(:) :: qq       ! Column specific humidity (kg/kg)
       real, intent(in), dimension(:) :: uu       ! Column zonal wind (m/s)
       real, intent(in), dimension(:) :: vv       ! Column meridional wind (m/s)
       real, intent(in), dimension(:) :: ww       ! Column vertical wind (m/s)
       real, intent(in), dimension(:) :: th       ! Column potential temperature (K)
       real, intent(in), dimension(:) :: pr       ! Column pressure (Pa)
    
       real, intent(out) :: UML        ! Mixed-layer horizontal velocity magnitude
       real, intent(out) :: Um         ! Mixed-layer zonal wind
       real, intent(out) :: Vm         ! Mixed-layer meridional wind
       real, intent(out) :: betaFC     ! Free convection beta parameter
       real, intent(out) :: DthFC      ! Free-convection plume excess potential temperature (K)
       real, intent(out) :: zFC        ! Free convection height above ground level (m)
       real, intent(out) :: pFC        ! Free convection pressure (Pa)
       real, intent(out) :: Bflux      ! Net buoyancy flux required for the plume to reach zFC (m^4.s^-3)
       real, intent(out) :: Hflux      ! Net heat flux required for the plume to reach zFC (J.s^-1 = W)
    
    
    ARGUMENTS
    ---------
    
    
    RETURN
    ------
    UML,Um,Vm,betaFC,DthFC,zFC,pFC,Bflux,Hflux
    """
    
    #Prcntg to pass an index referring to the cosphi values here:
    cosphi = [1.0, 0.8054, 0.6870, 0.5851, 0.4919, 0.4040, 0.3197, 0.2379, 0.1577, 0.0786, 0.0,
              -0.0786, -0.1577, -0.2379, -0.3197, -0.4040, -0.4919, -0.5851, -0.6870, -0.8054, -1.0]
    
    pd = len(TT) # how many pressure levels
    
    # Find the lowest pressure level above the land surface
    ksfc = numpy.sum(pr > psfc) - 1
    print("debug: ksfc",ksfc)
    #write(6,*) "ksfc",ksfc
    
    # Find the height index above the mixed-layer LCL (assumes pressure indices begin
    # with 1 at the lowest model level)
    #  Need to force the code to consider more than one pressure level so that a fog
    #  layer does not register as the LCL.  In the future this might need to be set
    #  to a minimum height (e.g., 250 m) when there are more frequent pressure
    #  levels.  Option added 13-09-2018.

    #write(6,*) "Calculating LCL:"
    #write(6,*) "Theta (C),      q (g/kg),     P_LCL (hPa),    T_LCL (C),  pr(k) (hPa), k" 
    #do k=ksfc+1,pd # skip first 2 levels in case of fog issue
    for k in range(ksfc+3,pd):
         recipk = 1.0/float(k)
         thmean = numpy.sum(th[:k])*recipk
         qmean  = numpy.sum(qq[:k])*recipk
         PLCL, TLCL = LCL(thmean,qmean,PLCL,TLCL)
         #write(6,*) thmean-273.15,qmean*1000,PLCL*0.01,TLCL-273.15,pr(k)*0.01,k
         if (PLCL > pr(k)):
             kLCL = k
             break
         #write(6,*) thmean-273.15,qmean*1000,PLCL*0.01,TLCL-273.15,pr(k)*0.01,k
    #write(6,*) "KLCL, pr(KLCL), TT(KLCL)",KLCL, pr(KLCL), TT(KLCL)
    
    
def LCL(th,qq):
    """
    ! Use Bolton's (1980) formula to calculate the lifting condensation level (LCL)
    ! temperature and pressure.
    
    real, intent(in) :: th     ! potential temperature
    real, intent(in) :: qq     ! specific humidity

    real, intent(out):: pLCL   ! pressure at the LCL
    real, intent(out):: TLCL   ! temperature at the LCL
    """
    
    #! Local variables
    #   real :: ev                 ! vapour pressure (Pa)
    #   real :: rr                 ! mixing ratio (kg/kg)

    # Calculate the vapour pressure valid at pressure = P0
    ev = P0/( (1-epsln) + epsln/qq)

    # Calculate the mixing ratio
    rr = qq/(1 - qq)

    # Calculate the LCL temperature.  Here we can use theta in place of the air
    # parcel temperature because the vapour pressure was evaluated at p = p0.
    # Note, the vapour pressure is multiplied by 0.01 to convert to units of hPa.
    TLCL = 2840.0/( 3.5*numpy.log(th) - numpy.log(ev*0.01) - 4.805 ) + 55.0

    # Calculate the LCL pressure
    pLCL = P0*(TLCL/th)**( Cpd/(Rd*(1.0-0.24*rr)) )

    return pLCL, TLCL