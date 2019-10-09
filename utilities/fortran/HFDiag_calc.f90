!
!  This program calculates the net heatflux required for
!  pyrocumulonimbus formation, from an atmospheric sounding.  
!  It will be called within a python script.  
!
!  For initial testing a single profile T, q, u, v and w will be 
!  specified in the main program body, as will all the parameters.
!
!  K. Tory Jan 2018.
!
!------------------------------------------------------

module physical_constants
   implicit none
   save

   real, parameter :: Cpd=1005.7           ! Specific heat of dry air (J/kg/K)
   real, parameter :: Rd=287.04            ! Gas constant for dry air (J/kg/K)
   real, parameter :: LV=2.501e6           ! Latent heat of vapourisation (J/kg)
   real, parameter :: P0=1.0e5             ! Standard pressure (Pa)
   real, parameter :: grav=9.8             ! Acceleration due to gravity (m/s^2)
   real, parameter :: epsln=0.622          ! Ratio of gas constants of dry air and water
   real, parameter :: pi=3.1415926536      ! 

end module physical_constants

module subroutines
   implicit none  
   save

contains
   subroutine heat_flux_calc(TT,qq,uu,vv,ww,th,pr,pd,zsfc,psfc,Tsfc, &
                             phi,DbSP,ni,nj,zp,beta_e,Pmin,betaSPmin,Wmin,Umin,Prcntg, &
                             UML,Um,Vm,betaFC,zFC,pFC,Bflux,Hflux)
   use physical_constants
   implicit none

! Declare input constants
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
 
! Declare input arrays
   real, intent(in), dimension(:) :: TT       ! Column temperature (K)
   real, intent(in), dimension(:) :: qq       ! Column specific humidity (kg/kg)
   real, intent(in), dimension(:) :: uu       ! Column zonal wind (m/s)
   real, intent(in), dimension(:) :: vv       ! Column meridional wind (m/s)
   real, intent(in), dimension(:) :: ww       ! Column vertical wind (m/s)
   real, intent(in), dimension(:) :: th       ! Column potential temperature (K)
   real, intent(in), dimension(:) :: pr       ! Column pressure (Pa)

! Declare output variables
   real, intent(out) :: UML        ! Mixed-layer horizontal velocity magnitude
   real, intent(out) :: Um         ! Mixed-layer zonal wind
   real, intent(out) :: Vm         ! Mixed-layer meridional wind
   real, intent(out) :: betaFC     ! Free convection beta parameter
   real, intent(out) :: zFC        ! Free convection height above ground level (m)
   real, intent(out) :: pFC        ! Free convection pressure (Pa)
   real, intent(out) :: Bflux      ! Net buoyancy flux required for the plume to reach zFC (m^4.s^-3)
   real, intent(out) :: Hflux      ! Net heat flux required for the plume to reach zFC (J.s^-1 = W)

! Declare local variables
   real :: recipk                  ! Reciprocal of the number of k levels used to sum fields
   real :: thmean                  ! Average value of theta (potential temperature, K)
   real :: qmean                   ! Average value of specific humidity (kg/kg)
   real :: PLCL                    ! Pressure at the lifting condensation level (Pa)
   real :: TLCL                    ! Temperature at the LCL (K)
   integer :: kLCL                 ! k-index immediately below the LCL
   integer :: ksfc                 ! k-index immediately below the surface
   real :: plower                  ! Pressure at the base of the height increment
   real :: Tlower                  ! Temperature at the base of the height increment


! Declare local miscellaneous variables
   real :: sum_wt                  ! Sum of the weighting terms used to normalise the weighting factors
   integer :: k                    ! Vertical profile (column) index
   integer :: i                    ! Saturation Point curve index
   integer :: j                    ! zb iteration index
   real :: zh                      ! Half height level used in calculating the weighting factors
   real :: thWML                   ! Weighted Mixed-layer mean potential temperature 
   real :: qqWML                   ! Weighted Mixed-layer mean specific humidity
   real :: WML                     ! Mixed-layer vertical velocity
   real :: s_star_max              ! Maximum saturation specific entropy of the column at or below Pmin
   integer :: imax                 ! SP curve index where the SP curve s* is first > s_star_max
   real :: alpha                   ! Linear interpolation factor (fraction of grid space where desired value resides)
   real :: betaSP                  ! beta value on SP curve corresponding to the position where s* = s_star_max
   real :: thFC                    ! Free convection potential temperature (cerresponding to betaFC)
   real :: qqFC                    ! Free convection specific humidity (corresponding to betaFC)
   real :: rhoFC                   ! Free convection density (corresponding to betaFC)
   real :: TFC                     ! Free convection temperature (corresponding to betaFC)
   integer :: kFC                  ! k-index of the pressure level just above pFC
   integer :: jconv                ! j-index when zb has converged
   real :: diff                    ! Difference between succesive iterations of zb
   integer :: kl                   ! k-index below the SP curve, temperature trace intersection
   integer :: kup                  ! k-index immediately above Pmin
   real :: AA,BB                   ! Constants for linear function of pressure with theta
   real :: thenv                   ! environment potential temperature
   real :: qenv                    ! environment specific humidity


! Declare local vertical profile arrays
   real, dimension(0:pd) :: zz     ! Height array at each pressure level with zsfc in the k=0 position
   real, dimension(pd)   :: wt     ! Weighting factor used to calculate a weighted mean th and qq 
                                   !   representing the relative fractions of th and qq entrained into a plume

! Declare local arrays along the SP curve
   real, dimension(ni) :: sSP      ! Saturation specifc entropy. (Column s* is temporarily stored in sSP(1))
   real, dimension(ni) :: bSP      ! Buoyancy parameter
   real, dimension(ni) :: pSP      ! Pressure
   real, dimension(ni) :: TSP      ! Temperature
   real, dimension(ni) :: thpl     ! Plume potential temperature 
   real, dimension(ni) :: qqpl     ! Plume specific humidity

! Declare local iteration arrays
   real, dimension(nj) :: zb       ! Effective plume centreline height (m)
   real, dimension(nj) :: f        ! Newton's method f(j)=0
   real, dimension(nj) ::fdash     ! Newton's Method first derivative of f(j) (f'(j) )

! Arrays for plume circle cuts
!   real, dimension(10) :: X        ! Fraction of plume circle above zFC
   real, dimension(21) ::cosphi    ! Fraction of plume radius between the plume circle origin and zFC
   real :: bcosphi                 ! beta_e*cosphi(Prcntg)

! Fill X and cosphi
!   data X/0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, &
!          0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.1/
!  Now using Prcntg to pass an index referring to the cosphi values here:
   data cosphi/1.0, 0.8054, 0.6870, 0.5851, 0.4919, 0.4040, 0.3197, 0.2379, 0.1577, 0.0786, 0.0, &
               -0.0786, -0.1577, -0.2379, -0.3197, -0.4040, -0.4919, -0.5851, -0.6870, -0.8054, -1.0/
!   write(6,*) "Inside subroutine heat_flux_calc"

! Find the lowest pressure level above the land surface
   ksfc = 0
   do k=1,pd
!     write(6,*) "pr(k),psfc,k",pr(k),psfc,k
     if(pr(k) <= psfc) then
        ksfc = k-1
        exit
     endif
   enddo
!   write(6,*) "ksfc",ksfc

! Find the height index above the mixed-layer LCL (assumes pressure indices begin
! with 1 at the lowest model level)
!  Need to force the code to consider more than one pressure level so that a fog
!  layer does not register as the LCL.  In the future this might need to be set
!  to a minimum height (e.g., 250 m) when there are more frequent pressure
!  levels.  Option added 13-09-2018.

   do k=ksfc+1,pd
     recipk = 1.0/real(k)
     thmean = sum(th(1:k))*recipk
     qmean  = sum(qq(1:k))*recipk
     if (k > ksfc+2) then  ! Added to deal with fog possibility
       call LCL(thmean,qmean,PLCL,TLCL)
!       write(6,*) thmean-273.15,qmean*1000,PLCL*0.01,TLCL-273.15,pr(k)*0.01,k
       if (PLCL > pr(k)) then
         kLCL = k
         exit
       endif
     endif
   enddo

!   write(6,*) "thmean,qmean,PLCL,TLCL,kLCL",thmean-273.15,qmean*1000,PLCL*0.01,TLCL-273.15,kLCL

! Calculate the weighting factors for a weighted mixed-layer LCL
   ! First fill the height array zz
   zz(0) = zsfc
   plower = psfc
   Tlower = Tsfc
   do k=1,pd
     if(pr(k) > psfc) then
       zz(k) = zsfc
     else
       zz(k) = zz(k-1) + log( (plower)/pr(k) )*Rd*(TT(k)+Tlower)/(2*grav)
       plower = pr(k)  ! Reset for next iteration
       Tlower = TT(k)
     endif
   enddo
!   write(6,*) "zz",zz
   ! Now all sub-ground height levels are set to the topographic height
   ! and ksfc is the k-index of the pressure level just below the sfc
   ! 
   wt = 0.0
   sum_wt = 0.0
   do k=ksfc+1,kLCL
     zh = (zz(k)+zz(k-1))*0.5
     if (k /= 1) then
       wt(k) = (zh**zp)*(pr(k-1)-pr(k)) 
     else
       wt(k) = (zh**zp)*(psfc-pr(k))
     endif
     sum_wt = sum_wt + wt(k)
!     write(6,*) "zh,wt(k)",zh,wt(k),k
   enddo 
   ! Express the weighting factor as a fraction of the total wt
   wt = wt/sum_wt
!   write(6,*) "wt",wt
   ! Here the weighting factor wt(k) represents the layer between k and
   ! k-1

! Calculate the weighted mixed-layer potential temperature and specific
! humidity. 
!  These calculations ignore the surface/near-surface values.
   thWML = 0.0
   qqWML = 0.0
   do k = ksfc+1,kLCL
     if (k /= 1) then
       thWML = thWML + wt(k)*(th(k)+th(k-1))*0.5 
       qqWML = qqWML + wt(k)*(qq(k)+qq(k-1))*0.5
     else
       thWML = th(k)*wt(k)
       qqWML = qq(k)*wt(k)
     endif
!     write(6,*) th(k),qq(k),k
   enddo
!   write(6,*) "thWML,qqWML",thWML-273.15,qqWML*1000

! Calculate the  mixed-layer horizontal (UML) and vertical (WML) winds.
! Here UML is the magnitude of the average vector wind (not the average wind
! speed).  Thus components on different levels can cancel.
!  These calculations ignore the surface/near-surface values.
   recipk = 1.0/real(kLCL-ksfc)
   WML = sum(ww(ksfc+1:kLCL))*recipk
   Um = sum(uu(ksfc+1:kLCL))*recipk
   Vm = sum(vv(ksfc+1:kLCL))*recipk
   UML = sqrt(Um**2 + Vm**2)
   UML = amax1(UML,Umin)   
!   write(6,*) "UML,WML",UML,WML

! Determine the maximum value of saturation specific entropy (s*) 
! on the pressure levels above the LCL up to the specified upper 
! pressure level(Pmin).  
! 
!---------------------------------------------------------------------
!  14-08-2018
! Dean Sgarbossa suggested using the pressure level where the 
! plume parcel temperature cools to -20 C (an electrification limit)
! to define the upper pressure level.
!
! Using the pressure level at which the environment T first
! drops below -20 C should work.  Need to pass a variable called
! Telec with value 253.15 K (electrification temperature limit). 
!
!  28-03-2019
! An alternative is to calculate the pressure at which the environment
! temperature is -20 C and pass this value as Pmin. Adjustments to the code have
! been added below so that Pmin is tested for maximum saturation specific
! entropy.
!---------------------------------------------------------------------
!
! Begin by finding the intersection of the SP curve with the 
! environment temperature trace.  s* at this point should be the 
! minimum possible s*_max.  To do this we need to identify the
! pressure levels above and below the intersection, and fill an 
! array of P and theta at very small increments between between
! these two levels.  By progressing up the SP curve increment by
! increment until both Pinc > PSP, and th_inc < theta_pl, then we
! have arrived at the first increment slightly beyond the intersection.
! Here the first s* should be calculated, and the increment index 
! recorded, since the looping along the SP curve below can start 
! from here.

   call find_k_index_below_intersection(th,qq,TT,pr,pd,thWML,qqWML,phi,kl)

!   write(6,*) "kl,pr(kl),TT(kl)",kl,pr(kl)*0.01,TT(kl)-273.15

! Calculate linear function thenv = AA*pr + BB
   AA = (th(kl+1) - th(kl))/(pr(kl+1) - pr(kl))
   BB = th(kl) - AA*pr(kl)

! Proceed along the SP curve in small increments.  Determine the pressure
! at each step, and find the environment theta at that pressure level. 
! When the theta on the SP curve first exceeds the environment, the 
! intersecton has been exceeded.  Assume this is close enough to the 
! true intersection position.  Store the SP curve increment in 'imax'.
! Calculate s* at this position and set it to s_star_max.
   do i = 1,ni
     bSP(i) = (i-1)*DbSP
     thpl(i) = (1 + bSP(i))*thWML
     qqpl(i) = qqWML + bSP(i)*phi*thWML
     call LCL(thpl(i),qqpl(i),pSP(i),TSP(i))
     thenv = AA*pSP(i) + BB 
!     write(6,*) "thpl,pSP,thenv",thpl(i)-273.15,pSP(i)*0.01,thenv-273.15
     if(thpl(i) > thenv) then
       imax = i
       call sat_spec_entr(thpl(i),qqpl(i),pSP(i),s_star_max)
       exit
     endif
     if(i == ni) then
       write(6,*) "WARNING: variable 'ni' is too small: ni =",ni
       write(6,*) "Increase the value in the driving script."
       write(6,*) "Program will stop!"
       stop
     endif
   enddo

!   write(6,*) "s_star_max",s_star_max

! Calculate s* at all the levels above kl, and search for the maximum
! value below Pmin.
! Store s* in the first position of the sSP array. 
   do k = kl+1,pd
     if(pr(k) > Pmin) then
       call sat_spec_entr(th(k),qq(k),pr(k),sSP(1))
!       write(6,*) sSP(1),pr(k)*0.01,k
       if(sSP(1) > s_star_max) s_star_max = sSP(1)
     else
       kup = k
       exit
     endif
   enddo
!   write(6,*) "Old s_star_max = ",s_star_max

! Calculate s* at Pmin and test for s_star_max (Added 2019-03-28)
   ! Calculate theta and q at Pmin and store in thenv and qenv
   alpha = (Pmin - pr(kup))/(pr(kup-1) - pr(kup))
   thenv = alpha*( th(kup-1) - th(kup) ) + th(kup)
   qenv  = alpha*( qq(kup-1) - qq(kup) ) + qq(kup)
   call sat_spec_entr(thenv,qenv,Pmin,sSP(1))
   if(sSP(1) > s_star_max) s_star_max = sSP(1)
!   write(6,*) "New s_star_max = ",s_star_max
!   write(6,*) "Pmin,pr(kup),pr(kup-1),alpha",Pmin,pr(kup),pr(kup-1),alpha
!   write(6,*) "thenv,th(kup),th(kup-1)",thenv,th(kup),th(kup-1)
!   write(6,*) "qhenv,qq(kup),qq(kup-1)",qenv,qq(kup),qq(kup-1)

!   write(6,*) "s_star_max",s_star_max

! Loop through increments along the SP curve starting at the WML-LCL

   sSP = 0.0   ! Reset sSP array
   do i = imax,ni
     bSP(i) = (i-1)*DbSP
     thpl(i) = (1 + bSP(i))*thWML
     qqpl(i) = qqWML + bSP(i)*phi*thWML
     call LCL(thpl(i),qqpl(i),pSP(i),TSP(i))
     call sat_spec_entr(thpl(i),qqpl(i),pSP(i),sSP(i))
!     write(6,*) sSP(i),thpl(i)-273.15,qqpl(i)*1000,pSP(i)*.01,bSP(i),i
     if (sSP(i) >= s_star_max) then   ! Maximum i-index found
       imax = i
       exit
     endif
   enddo
!   write(6,*) "sSP,bSP,pSP,thpl,qqpl",sSP(imax),bSP(imax),pSP(imax)*0.01,thpl(imax)-273.15,qqpl(imax)*1000,imax

! Use linear interpolation to find specific values of betaSP, PrSP that
! more closely match where sSP = s_star_max
   if(imax > 1) then
     alpha = (s_star_max - sSP(imax-1))/(sSP(imax) - sSP(imax-1))
     betaSP = bSP(imax-1) + alpha*(bSP(imax) - bSP(imax-1))
   else
     betaSP = bSP(1)
   endif

! Add the minimum buoyancy to betaSP to get betaFC (free convection
! beta).  This is to account for some dilution of the rising cloud.  It
! could be defined by how dry the air is above the LCL (a future
! project perhaps).
   betaFC = betaSP + betaSPmin

! Calculate free convection pressure (pFC) and density (rhoFC)
   thFC = (betaFC+1)*thWML
   qqFC = qqWML + betaFC*phi*thWML
   call LCL(thFC,qqFC,pFC,TFC)
   rhoFC = pFC/(Rd*TFC)

! Use linear interpolation to get the free convection height (zFC)
! First find the pressure level above the free convection level
! Start the search below kLCL to account for cases of bFC < 0
! Subtract zfc to get height above ground level (added 29-01-2019)
   do k = ksfc+1,pd
     if(pr(k) < pFC) then  ! pressure level above pFC found
       kFC = k
       exit
     endif
   enddo
!   write(6,*) "kFC,thFC,qqFC,rhoFC",kFC,thFC,qqFC,rhoFC
   alpha = (pFC - pr(kFC-1))/(pr(kFC) - pr(kFC-1))
   zFC = zz(kFC-1) + alpha*(zz(kFC) - zz(kFC-1)) - zsfc

! Begin the net Buoyancy flux calculation (Bflux = pi*BB)
! Assume Briggs model for now.  May consider a hybrid Briggs/MTT model
! in the future.  If the mean mixed-layer vertical motion (WML) is zero
! the calculation is simple and direct.  If WML is non-zero an iterative
! method is needed to solve for zb, an adjusted height that takes into account
! the difference in height between the plume centreline and the plume top and
! the impact of vertical motion on the plume.

   bcosphi = beta_e*cosphi(Prcntg)
   zb(1) = zFC/(1+bcosphi)
   if (sqrt(WML**2) > Wmin) then ! W is significant
     jconv = 0
     do j=1,nj
       f(j) = zFC - (1+bcosphi)*zb(j) - sqrt( (2.0*zb(j))/(3.0*grav*betaFC) )*WML
       fdash(j) = -(1+bcosphi) - 0.5*sqrt( (2.0)/(3.0*grav*betaFC*zb(j)) )*WML
       zb(j+1) = zb(j) - f(j)/fdash(j)
       diff = sqrt( ( (zb(j+1)-zb(j))/zb(j) )**2  )
!       write(6,*) "diff",diff,j
       if (diff < 0.002) then ! Convergence achieved. Place converged value in zb(1) 
         zb(1) = zb(j+1)
         jconv = j 
         exit
       endif
     enddo
     if (jconv == 0) then
       write(6,*) "WARNING: zb did not converge after ",j,"iterations"
       stop
     endif
   endif

! Calculate BB, Bflux and Hflux
! It's not clear what value of rho is most appropriate.  Is it the density of
! the combustion gases, or the plume at some elevated level?  This might depend
! on whether the heat flux or buoyancy flux is conserved in the plume.  It may
! be that the density does not change much (pressure and temperature decrease at
! a similar rate) in which case both BB and HH are conserved, and we can choose
! any value of rho that we have suitable data for.  
   Bflux = pi*grav*betaFC*UML*(beta_e*zb(1))**2
   Hflux = (rhoFC*Cpd*thWML/grav)*Bflux
 
!   write(6,*) "UML, WML,Umin,Wmin",UML,WML,Umin,Wmin

!   write(6,*) "Buoyancy Froude Number =",UML/sqrt(betaFC*grav*zb(1))
!   write(6,*) "dzb/dx =",sqrt(2./3.*betaFC*grav*zb(1))/UML

   end subroutine heat_flux_calc

!----------------------------------------------------------

   subroutine find_k_index_below_intersection(th,qq,TT,pr,pd,thWML,qqWML,phi,kl)
   implicit none

   integer, intent(in) :: pd       ! Number of elements in the pressure (vertical) dimension
   real, intent(in), dimension(pd) :: TT       ! Column temperature (K)
   real, intent(in), dimension(pd) :: qq       ! Column specific humidity (kg/kg)
   real, intent(in), dimension(pd) :: th       ! Column potential temperature (K)
   real, intent(in), dimension(pd) :: pr       ! Column pressure (Pa)

   real, intent(in) :: thWML     ! Weighted mixed-layer potential temperature (K)
   real, intent(in) :: qqWML     ! Weighted mixed-layer specific humidity (kg/kg)
   real, intent(in) :: phi       ! Fire moisture to heat ratio (kg/kg/K)
   integer, intent(out) :: kl    ! k-index below the SP curve, temperature trace intersection

! Local arrays
   real, dimension(pd) :: beta   ! Value of beta parameter where the model level theta surface 
                                    ! intersects the SP curve
   real :: qbeta                 ! Specific humidity on the SP curve at beta
   real :: Pbeta                 ! Pressure on the SP curve at beta
   real :: Tbeta                 ! Temperature on the SP curve at beta
   
   integer :: k                                ! Loop index

! Fill beta array
   do k=1,pd
     beta(k) = th(k)/thWML - 1
!     write(6,*) "beta",beta(k),k
   enddo

! Find the first guess lower k-index (kl)
   do k=1,pd
     if(beta(k) > 0) then
       kl = k
       exit
     endif
   enddo
!   write(6,*) "kl",kl

! Find the first lower k-index where Tbeta is greater than the temperature on the upper k-index.  
   do k=kl,pd
     qbeta = qqWML + beta(k+1)*phi*thWML
     call LCL(th(k+1),qbeta,Pbeta,Tbeta)
!     write(6,*) "Tbeta,Pbeta,TT(k+1)",Tbeta-273.15,Pbeta*0.01,TT(k+1)-273.15,k+1
     if(Tbeta > TT(k+1)) then
       kl = k
       exit
     endif
   enddo

   end subroutine find_k_index_below_intersection

!----------------------------------------------------------


   subroutine LCL(th,qq,pLCL,TLCL)
   use physical_constants
   implicit none

! Use Bolton's (1980) formula to calculate the lifting condensation level (LCL)
! temperature and pressure.

! Array declaration
   real, intent(in) :: th     ! potential temperature
   real, intent(in) :: qq     ! specific humidity
   real, intent(out):: pLCL   ! pressure at the LCL
   real, intent(out):: TLCL   ! temperature at the LCL

! Local variables
   real :: ev                 ! vapour pressure (Pa)
   real :: rr                 ! mixing ratio (kg/kg)

! Calculate the vapour pressure valid at pressure = P0
   ev = P0/( (1-epsln) + epsln/qq)

! Calculate the mixing ratio
   rr = qq/(1 - qq)

! Calculate the LCL temperature.  Here we can use theta in place of the air
! parcel temperature because the vapour pressure was evaluated at p = p0.
! Note, the vapour pressure is multiplied by 0.01 to convert to units of hPa.
   TLCL = 2840.0/( 3.5*log(th) - log(ev*0.01) - 4.805 ) + 55.0

! Calculate the LCL pressure
   pLCL = P0*(TLCL/th)**( Cpd/(Rd*(1.0-0.24*rr)) )

   end subroutine LCL

!----------------------------------------------------------
   subroutine sat_spec_entr(th,qq,pr,s_star)
   use physical_constants
   implicit none

! Array declaration
   real, intent(in) :: th     ! potential temperature
   real, intent(in) :: qq     ! specific humidity
   real, intent(in) :: pr     ! pressure
   real, intent(out):: s_star ! saturation specific entropy

! Local variables
   real :: e_star             ! saturation vapour pressure (Pa)
   real :: pd_star            ! saturation partial pressure of dry air (Pa)
   real :: r_star             ! saturation mixing ratio (kg/kg)
   real :: TT                 ! temperature (K)

! Calculate mixing ratio and store in r_star temporarily
   r_star = qq/(1 - qq)

! Calculate e_star
   TT = th*(pr/P0)**(Rd/Cpd*(1 - 0.24*r_star))
   e_star = 611.2*exp( (17.67*TT- 4826.56)/(TT - 29.65) )
   
! Calculate pd_star
   pd_star = pr - e_star

! Calculate r_star
   r_star = epsln*e_star/pd_star

! Calculate s_star
   s_star = Cpd*log(TT) + Lv*r_star/TT - Rd*log(pd_star)

   end subroutine sat_spec_entr


   


end module subroutines

!------------------------------------------------------
program Heatflux
   use subroutines
   use physical_constants
   implicit none

! Declare input constants
   integer, parameter  :: pd=22       ! Number of elements in the pressure (vertical) dimension
   real :: zsfc        ! Surface height (m)
   real :: psfc        ! Surface pressure (Pa)
   real :: Tsfc        ! Surface temperature (K)
   real :: phi         ! Ratio of fire moisture to heat (
   real :: DbSP        ! Beta increment along the saturation point curve
   integer :: ni       ! Number of increments along the SP curve
   integer :: nj       ! Number of iterations for solving plume centre-line height
   real :: beta_e      ! Briggs entrainment paramter (Briggs uses 0.6, we may use up to 0.9)
   real :: zp          ! The power to which zz is raised in calculating wt. (2/3 for MTT, 1.0 for Briggs)
   real :: Pmin        ! The minimum pressure the moist plume needs to rise to for 
                                   !   pyroCb to be considered.  Probably choose
                                   !   300 to 250 hPa, 
                                   !   just below the tropopause.
   real :: betaSPmin   ! A minimum value of beta to be added to betaSP, to account 
                                   !   for buoyancy losses from entrainment etc.
                                   !   above the condensation level
   real :: Wmin        ! Minimum value of mixed-layer wind speed to be considered to influence zb
   real :: Umin        ! Minimum value of mixed-layer horizontal wind speed
   integer :: Prcntg   ! Indec for specifying percentage of plume required to reach zFC


! Declare input arrays
   real, dimension(pd) :: TT       ! Column temperature (K)
   real, dimension(pd) :: qq       ! Column specific humidity (kg/kg)
   real, dimension(pd) :: RH       ! Column relative humidity
   real, dimension(pd) :: uu       ! Column zonal wind (m/s)
   real, dimension(pd) :: vv       ! Column meridional wind (m/s)
   real, dimension(pd) :: ww       ! Column vertical wind (m/s)
   real, dimension(pd) :: th       ! Column potential temperature (K)
   real, dimension(pd) :: pr       ! Column pressure (Pa)

! Declare output variables
   real :: betaFC     ! Free convection beta parameter
   real :: zFC        ! Free convection height (m)
   real :: pFC        ! Free convection pressure (Pa)
   real :: Bflux      ! Net buoyancy flux required for the plume to reach zFC (m^4.s^-3)
   real :: Hflux      ! Net heat flux required for the plume to reach zFC (J.s^-1 = W)
   real :: UML        ! Mixed-layer horizontal velocity magnitude
   real :: Um         ! Mixed-layer zonal wind
   real :: Vm         ! Mixed-layer meridional wind

! Local variables
   real :: e_star     ! saturation vapour pressure (Pa)
   real :: ee         ! vapour pressure (Pa)
   real :: rr         ! mixing ratio (kg/kg)
   integer :: k       ! Loop variable

! Fill arrays

   data pr/1000, 975, 950, 925, 900, 850, 800, 750, 700, 600, 500, 450, 400, 350, &
           300, 275, 250, 225, 200, 175, 150, 100/

!! Sir Ivan 20170212 0000
!   data TT/313.1292, 311.6246, 310.0886, 308.4834, 306.1262, 301.5018, 296.8003, &
!           291.6119, 286.3469, 275.4677, 269.8451, 265.1353, 258.8537, 251.1692, &
!           243.1011, 238.6128, 233.7222, 227.6593, 223.1471, 217.1296, 209.8734, &
!           195.693/
!   data RH/10.07302, 10.07302, 10.07302, 12.52325, 13.88663, 17.22936, 22.49171, &
!           31.03843, 43.81114, 45.65998, 1.251922, 1.029724, 8.305652, 8.234468, &
!           3.300921, 17.794, 21.5922, 23.79617, 5.658082, 9.43286, 4.901436, &
!           15.98808/
!   data uu/7.664062, 7.664062, 7.664062, 13.63928, 14.98028, 17.5642, 21.29099, &
!           21.0274, 16.52952, 12.69814, 19.04091, 20.75769, 22.05462, 23.95693, &
!           25.82995, 24.06316, 23.99809, 25.26415, 23.13675, 27.3345, 27.72713, &
!           17.034/
!   data vv/-1.656248, -1.656248, -1.656248, -2.65053, -2.589676, -2.294632, -2.797358, &
!           -4.943251, -8.957104, -12.75691, -7.239407, -7.144512, -6.573563, -5.138931, &
!           -9.482964, -13.82598, -15.87357, -11.12967, -9.575665, -13.81601, -10.48209, &
!           -9.852392/
!   data ww/0.003222659, 0.003222674, 0.003222644, 0.04492188, 0.04034775, 0.03677505, &
!           0.07510799, 0.08570462, 0.05466044, 0.03466278, -0.01228929, -0.03957504, &
!          -0.06792116, -0.04204297, -0.01746798, -0.03679788, -0.01223761, 0.02117497, &
!           0.02292418, 0.00498575, -0.01135153, 0.006939173/

! Set surface variables Sir Ivan 20170212 0000
   zsfc = 492.625
   psfc = 94657.5
   Tsfc = 318.0625

! Sir Ivan 20170212 1200
   data TT/297.4067, 295.9777, 294.5247, 292.9023, 290.7185, 293.0745, 292.9617, &
           289.9343, 284.7459, 274.384, 265.941, 261.2345, 256.1496, 248.879, &
           242.9845, 239.0938, 234.1526, 229.0922, 223.8137, 217.1339, 210.5804, &
           199.3855/
   data RH/37.88145, 37.88145, 37.88406, 41.29988, 42.69291, 33.30973, 41.16108, &
           34.16614, 43.67365, 74.46993, 62.73916, 46.99235, 65.02315, 64.0191, &
           34.54897, 26.02923, 12.38063, 1.941115, 1.265086, 6.369132, 5.15822, &
           9.187106/
   data uu/2.195314, 2.195312, 2.195312, 4.154184, 4.88792, 13.41563, 17.93077, &
           18.73009, 20.2985, 21.04253, 20.80694, 25.43085, 29.38869, 29.83351, &
           32.70073, 32.48066, 33.18607, 35.63335, 36.50887, 39.11257, 40.0449, &
           23.30711/
   data vv/4.671875, 4.671877, 4.671875, 9.290928, 10.37664, 10.46318, 4.95372, &
           4.675549, 3.896252, -3.367939, -3.393467, -5.94738, -7.163158, -10.94949, &
          -11.1648, -13.52404, -13.28518, -10.66497, -12.24237, -14.87078, -14.84602, &
          -4.965075/
   data ww/0.001855478, 0.001855448, 0.00627172, 0.01837605, -0.001328826, -0.04578185, &
           -0.04078722, -0.02353448, -0.03373182, -0.06114858, 0.03424823, 0.1030626, &
            0.06273985, 0.06807274, 0.0008055568, -0.04409993, -0.04419011, -0.03981435,&
           -0.08006632, -0.1107656, -0.1095305, -0.0414575/

! Set surface variables Sir Ivan 20170212 1200
   zsfc = 492.625
   psfc = 95066.12
   Tsfc = 294.8594

!!  Yarloop 20160105 1200 <-- This was two days before the Yarloop event (Wrong day)
!   data TT/298.0977, 296.6047, 298.1581, 299.1079, 296.9154, 292.424, 288.242, &
!           283.9648, 279.7549, 272.2592, 262.9808, 257.6738, 250.5752, 242.6264, &
!           233.9513, 229.5312, 225.1123, 221.1917, 217.1844, 213.0707, 211.5655, &
!           207.9612/
!   data RH/60.83675, 62.66201, 50.88776, 42.56545, 45.9001, 53.69991, 66.02894, &
!           80.17151, 77.8252, 35.38151, 9.948997, 4.518318, 6.171204, 9.141273, &
!           14.21976, 17.2197, 17.67703, 12.01009, 15.30079, 17.89166, 12.133, &
!           4.237309/
!   data uu/0.6328125, 0.6574974, 0.003299713, -5.076874, -7.503613, -7.631676, &
!           -5.149014, -3.861458, -1.915451, 1.598793, 2.939644, 3.249092, &
!            4.645309, 6.532326, 9.419136, 9.670631, 11.84599, 16.21149, &
!            16.71908, 14.77872, 14.9489, 16.83636/
!   data vv/-1.726562, -2.871128, 3.761543, 7.096272, 7.150887, 6.528511, 5.324005, &
!            1.601677, -5.171127, 0.4387703, 4.386566, 1.174873, 1.01004, 1.313408, &
!            1.11692, 0.1325912, -1.316307, -1.461037, -2.157677, -5.439018, -7.384678, &
!            -4.187782/
!   data ww/0.0009765625, 0.01744103, 0.01086307, -0.01044703, -0.01159507, 0.005419493, &
!           -0.004862309, -0.03331232, -0.06031585, 0.02295333, 0.002535582, -0.03640324, &
!           -0.07250828, -0.07427537, -0.03842074, -0.002075315, 0.02259678, 0.02305108, &
!            0.04085571, 0.03228259, 0.002595007, -0.02993613/
!
!! Set surface variables   Yarloop 20160105 1200
!   zsfc = 140.875
!   psfc = 99453.38
!   Tsfc = 297.2344
!
!! Waroona 20160106 1200
!   data TT/306.0181, 304.5318, 302.7144, 300.987, 298.9988, 294.6937, 290.1977, 285.5171, &
!           280.421, 271.03, 262.2682, 257.4754, 250.2595, 242.7007, 234.4683, 229.7743, &
!           224.3409, 218.6008, 214.1742, 213.7108, 212.5852, 208.6684/
!   data RH/30.73079, 30.71032, 32.59037, 34.33974, 37.19561, 43.86477, 51.06636, 60.2851, &
!           72.31855, 34.39447, 4.528992, 3.12706, 6.540573, 9.326401, 11.15578, 14.42374, &
!           25.78455, 36.46635, 38.98806, 18.31938, 9.375191, 3.926466/
!   data uu/-5.195312, -12.20892, -16.12396, -16.15635, -14.80178, -12.59809, -10.24452,&
!           -7.889229, -5.880882, 2.645958, 1.904366, 1.189323, -0.5143204, -2.308258, &
!           -0.5762711, 1.140938, 1.563957, 1.726265, 5.323776, 6.041206, 5.047874, &
!            5.311501/
!   data vv/2.820312, 5.323475, 4.894508, 2.770092, 1.126724, -1.616142, -3.855846, &
!           -6.269821, -8.391594, -7.833473, -2.984509, -2.232765, -0.8088608, -0.8543167, &
!            0.3632584, 1.089752, 1.246864, 0.7059555, -0.04255676, -0.4358063, -0.09459686, &
!            2.208679/
!   data ww/-0.006445318, -0.1348806, -0.1738696, -0.1827778, -0.1728321, -0.1401746, &
!           -0.1089587, -0.04350758, 0.005830765, -0.004502356, -0.01326931, 0.01093227, &
!            0.03391498, 0.02372831, 0.02845699, 0.03598511, 0.02727336, 0.005504012, &
!           -0.01861954, -0.0136717, -0.00960511, 0.01235026/
! 
!! Set surface variables Waroona 20160106 1200
!   zsfc = 130.875
!   psfc = 99808.25
!   Tsfc = 303.9844

! Convert pressure to pascals
   pr = pr*100.0
!   ww = ww*1000.0  ! Artificially pump up ww to see what happens
!   uu = uu*2.0     ! Artificially pump up the background wind to see what happens
!   vv = vv*2.0

! Set constants
   phi = 6.0e-5  
   DbSP = 0.001
   ni = 101
   nj = 20
   beta_e = 0.4
   zp = 1.0
   Pmin = 35000 
   Wmin = 0.05
   Umin = 1.0
   betaSPmin = 0.002
   Prcntg = 7   ! 1 = 0.0 %, 2 = 5%, 3 = 10% and so on until 21 = 100%

!  Calculate specific humidity
   do k=1,pd
     ! saturation vapour pressure
     e_star = 611.2*exp( (17.67*TT(k)- 4826.56)/(TT(k) - 29.65) )
     ! vapour pressure
     ee = RH(k)*e_star/100.0
     ! specific humidity
     qq(k) = (epsln*ee)/( pr(k) - ee*(1.0-epsln) )
   enddo

!  Calculate potential temperature
   do k=1,pd
     rr = qq(k)/(1.0 - qq(k))
     th(k) = TT(k)*(p0/pr(k))**(Rd/Cpd*(1 - 0.24*rr))
!     write(6,*) rr,p0,(Rd/Cpd*(1 - 0.24*rr)),pr(k),th(k),k
!     write(6,*) pr(k)*0.01,TT(k)-273.15,qq(k)*1000.0
   enddo

!  write(6,*) TT(1),qq(1),uu(1),vv(1),ww(1),th(1),pr(1)
!  write(6,*) pd,zsfc,psfc,Tsfc
!  write(6,*) phi,DbSP,ni,nj,zp,beta_e,Pmin,betaSPmin,Wmin

  call heat_flux_calc(TT,qq,uu,vv,ww,th,pr,pd,zsfc,psfc,Tsfc, &
                             phi,DbSP,ni,nj,zp,beta_e,Pmin,betaSPmin,Wmin,Umin,Prcntg, &
                             UML,Um,Vm,betaFC,zFC,pFC,Bflux,Hflux)

!   write(6,*) "betaFC,zFC,pFC,Bflux,Hflux",betaFC,zFC,pFC,Bflux,Hflux

end program Heatflux
