module mod_thermo

real, parameter :: P_Rd  = 287.04
real, parameter :: P_Rv  = 461.50
real, parameter :: P_eps = P_Rd/P_Rv
real, parameter :: P_Cpd = 1004.64
real, parameter :: P_Cpv = 1870.0
real, parameter :: P_Lv  = 2.5e6

! Various thermodynamic functions.

contains

!** e_sat *******************************************************************
! Saturation vapour pressure in Pa at temperature T (K) and pressure p (Pa)
!
subroutine e_sat(T, p, esat, desdT)
   implicit none
   real, intent(in )           :: T, p
   real, intent(out)           :: esat 
   real, intent(out), optional :: desdT

   esat = (1.0007 + 3.46e-8*p)*611.21*exp(17.502*(T-273.15)/(T-32.18))
   if (present(desdT)) then
      desdT = 4217.46 * esat / (T - 32.18)**2
   endif
end subroutine e_sat
!
!** qsat ********************************************************************
! Given temperature T (K) and pressure p (Pa) returns the saturation specific
! humidity in kg/kg over pure water. Ref: Buck, JAM 20, 1981, 1527-1532.

      real function qsat(T, p)

      implicit none

      real p, T

      real esat

!     Saturation water vapour partial pressure in moist air over pure water

      call e_sat(T,p,esat)

!     Get saturation specific humidity

      qsat = P_eps * esat / (p - (1.0-P_eps)*esat)
!
      return
      end function qsat
!
!** qsatsea ******************************************************************
! Given temperature T (K) and pressure p (Pa) returns the saturation specific
! humidity in kg/kg over sea water.

      real function qsatsea(T, p)

      implicit none

      real p, T
      real esat

!     Saturation vapour pressure over sea water

      call e_sat(T,p,esat)
      esat = 0.98 * esat

!     Get saturation specific humidity

      qsatsea = P_eps * esat / (p - (1.0-P_eps)*esat)

      return
      end function qsatsea
!
!* r_sat ********************************************************************
! Saturation mixing ratio in kg/kg at temperature T (K) and pressure p (Pa)

subroutine r_sat(T, p, rs, drsdT)
   implicit none
   real, intent(in )           :: T, p
   real, intent(out)           :: rs 
   real, intent(out), optional :: drsdT

   real es, desdT

   if (.not. present(drsdT)) then
      call e_sat(T, p, es)
      rs = P_eps * es / ( p - es )
   else
      call e_sat(T, p, es, desdT)
      rs = P_eps * es / ( p - es )
      drsdT = P_eps * p * desdT / (p - es)**2
   endif
end subroutine r_sat

!* r_satsea *****************************************************************
! Saturation mixing ratio over sea water in kg/kg at temperature T (K) and 
! pressure p (Pa)
!
      real function r_satsea(T, p)
      implicit none
      real T, p
!
      real esat
!
      call e_sat(T, p, esat)
      esat = 0.98 * esat
      r_satsea = P_eps* esat / ( p - esat )
      return
      end function r_satsea

!* theta_l ************************************************************
subroutine theta_l(th, rtot, p, th_l, rsat, dthldth)
! Calculate the liquid water potential temperature, given T and rtot.
! Emanuel (1994) eq 4.5.15.

   implicit none
   real, intent(in)            :: th, rtot, p
   real, intent(out)           :: th_l
   real, intent(out), optional :: rsat, dthldth

   real chi, gamma, exner, T, rl, drsdT, drldth
   real a, b, c, dadth, dbdth, dcdth

   chi = (P_Rd + rtot*P_Rv) / (P_Cpd + rtot*P_Cpv)
   gamma = (rtot*P_Rv) / (P_Cpd + rtot*P_Cpv)

   exner = (p/1e5)**chi

   T = th * exner

   if (.not. present(dthldth)) then
       call r_sat(T, p, rsat)
       rl = max(rtot - rsat, 0.0)
    
       th_l = th * (1 - rl/(P_eps + rtot))**chi * (1 - rl/rtot)**(-gamma) &
           * exp(-P_Lv * rl / ((P_Cpd + P_Cpv*rtot)*T))
   else
       call r_sat(T, p, rsat, drsdT)
       if (rtot > rsat) then
          rl = rtot - rsat
          drldth = -drsdT*exner
       else
          rl = 0.0
          drldth = 0.0
       endif
       
       a = (1 - rl/(P_eps + rtot))**chi
       b = (1 - rl/rtot)**(-gamma)
       c = exp(-P_Lv * rl / ((P_Cpd + P_Cpv*rtot)*T))
    
       dadth = (-chi/(P_eps+rtot))*(1 - rl/(P_eps + rtot))**(chi-1)*drldth
       dbdth = (gamma/rtot)*(1 - rl/rtot)**(-gamma-1)*drldth
       dcdth = (P_Lv * (rl*exner - T*drldth) / ((P_Cpd + P_Cpv*rtot)*T**2)) * c
    
       th_l = th * a * b * c
       dthldth = a*b*c + th*(dadth*b*c + a*dbdth*c + a*b*dcdth)
   endif

end subroutine theta_l

!* inv_theta_l ********************************************************

subroutine inv_theta_l(thl, rtot, p, th, rsat)
! Inverts theta_l, using safe Newton-Raphson

   implicit none

   real, intent(in)  :: thl, rtot, p
   real, intent(out) :: th, rsat
!  real, intent(out), optional  rsat

   integer, parameter :: itmax = 40

   real chi, gamma, exner
   real err, dth, dthold
   real th0, th1, th2, T0, T1, T2, thl1, dthldth1
   real rsat0, rsat1, rsat2, rl0, rl1, rl2
   integer it

   chi = (P_Rd + rtot*P_Rv) / (P_Cpd + rtot*P_Cpv)
   gamma = (rtot*P_Rv) / (P_Cpd + rtot*P_Cpv)

   exner = (p/1e5)**chi

   ! Lower bound for th, guarantees thl0 <= thl
   th0 = thl
   T0 = th0 * exner
   call r_sat(T0, p, rsat0)
   rl0 = max(rtot - rsat0, 0.0)
   if (rl0 .eq. 0) then
      th = th0 
      rsat = rsat0
      return
   endif

   ! Upper bound for th, has thl <= thl2
   call r_sat(thl*exner, p, rsat2)
   rl2 = max(rtot - rsat2, 0.0)
   th2 = thl + (P_Lv/P_Cpd)*(rl2/exner)

   ! Initial guess for N-R part of algorithm
   th1 = 0.5*(th0 + th2)
   call theta_l(th1, rtot, p, thl1, rsat1, dthldth1) 

   dth = th2 - th1
   dthold = dth

   do it = 1, itmax
      ! Newton-Raphson step
      dthold = dth
      dth =  (thl1 - thl)/dthldth1
      th1 = th1 - dth
      ! Check within bounds and converging fast enough, else bisect
      if ( th1 < th0 .or. th1 > th2 .or. &
           abs(2*(thl1-thl)) > abs(dthold*(dthldth1 - 1)) ) then
         dth = 0.5*(th2 - th0)
         th1 = 0.5*(th0 + th2)
      endif
      call theta_l(th1, rtot, p, thl1, rsat1, dthldth1)
      err = abs(thl - thl1)
    
      if (err < 1e-4) then
         exit
      endif
    
      ! Update bounds
      if (thl1 < thl) then
         th0 = th1
      else
         th2 = th1
      endif
   enddo

   if (it < itmax) then
      th = th1
      rsat = rsat1
   else
      print *,'ERROR>>> Convergence failure in inv_theta_l'
      print *,'theta_l=',thl,' rtot=',rtot,' p=',p
      print *,th0,th1,th2
      stop
   endif

end subroutine inv_theta_l

end module mod_thermo

