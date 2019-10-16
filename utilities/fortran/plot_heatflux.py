# -*- coding: utf-8 -*-

'''
Extract and plot MARS/MSAS pressure level data (downloaded from /g/data/rr4)
'''

import matplotlib as mpl
mpl.use('Agg') # This goes before pyplot
import matplotlib.pyplot as plt # plotting resources
import matplotlib.colors as col
import matplotlib.ticker as tick
import numpy as np # for computing
np.seterr('ignore')
from netCDF4 import Dataset
import netCDF4 as ncd
import scipy.io as sio
import scipy.signal as signal
import datetime
from fortranbit import subroutines
#from matplotlib.backends.backend_pdf import PdfPages
import argparse
parser = argparse.ArgumentParser(description='Pass longitude and latitude plot window')
parser.add_argument('lonmin', metavar='LONMIN', type=float, nargs='?', default=110.0, help='Westernmost Longitude to be plotted')
parser.add_argument('lonmax', metavar='LONMAX', type=float, nargs='?', default=155.0, help='Easternmost Longitude to be plotted')
parser.add_argument('latmin', metavar='LATMIN', type=float, nargs='?', default=-45.0, help='Southernmost Latitude to be plotted')
parser.add_argument('latmax', metavar='LATMAX', type=float, nargs='?', default=-10.0, help='Northernmost Latitude to be plotted')

args = parser.parse_args()

lonmin=args.lonmin
lonmax=args.lonmax
latmin=args.latmin
latmax=args.latmax

print 'lonmin = ',lonmin
print 'lonmax = ',lonmax
print 'latmin = ',latmin
print 'latmax = ',latmax

datadir ='/home/toryk/WORK'
print(datadir)
figdir = datadir

# Map coordinates
map  = sio.loadmat('Australia_medium') # load matlab file
map_lon = map['Aust_medium'][:,0] 
map_lat = map['Aust_medium'][:,1]

ncfile = Dataset(datadir + '/fc_p.nc', 'r')
ncfile2 = Dataset(datadir + '/fc_sfc.nc', 'r')
ncfile3 = Dataset(datadir + '/an_sfc.nc', 'r')     # Added 09-09-2019 to read topography from the analysis file


print(ncfile)
print(ncfile2)
print(ncfile3)

#! Declare input constants
#   real :: phi         ! Ratio of fire moisture to heat (
#   real :: DbSP        ! Beta increment along the saturation point curve
#   integer :: ni       ! Number of increments along the SP curve
#   integer :: nj       ! Number of iterations for solving plume centre-line height
#   real :: beta_e      ! Briggs entrainment paramter (Briggs uses 0.4 for the internal plume)
#   real :: zp          ! The power to which zz is raised in calculating wt. (2/3 for MTT, 1.0 for Briggs)
#   real :: Pmin        ! The minimum pressure the moist plume needs to rise to for 
#                                   !   pyroCb to be considered.  Probably choose
#                                   !   300 to 250 hPa, 
#                                   !   just below the tropopause.
#   real :: betaSPmin   ! A minimum value of beta to be added to betaSP, to account 
#                                   !   for buoyancy losses from entrainment etc.
#                                   !   above the condensation level
#   real :: Wmin        ! Minimum value of mixed-layer wind speed to be considered to influence zb
#   real :: Umin        ! Minimum value of mixed-layer horizontal wind speed
#   integer :: Prcntg   ! Index for specifying percentage of plume required to reach zFC
#
#! Local variables
#   real :: e_star     ! saturation vapour pressure (Pa)
#   real :: ee         ! vapour pressure (Pa)
#   real :: rr         ! mixing ratio (kg/kg)

# Set constants
phi = 6.0e-5
DbSP = 0.001
ni = 1001
nj = 20
beta_e = 0.4
zp = 1.0
Pmin = 35000.0         # Minimum pressure the plume must rise to (Pa)
Tmin = 253.15         # Minimum temperature level the plume must rise to (K).
PorT = 'T'           # Set to 'T' if Tmin is to be used, otherwise set to 'P'
Wmin = 0.05
Umin = 1.0
betaSPmin = 0.002
Cpd=1005.7           # Specific heat of dry air (J/kg/K)
Rd=287.04            # Gas constant for dry air (J/kg/K)
LV=2.501e6           # Latent heat of vapourisation (J/kg)
p0=1.0e5             # Standard pressure (Pa)
grav=9.8             # Acceleration due to gravity (m/s^2)
epsln=0.622          # Ratio of gas constants of dry air and water
pi=3.1415926536
Prcntg = 7   # 1 = 0.0 %, 2 = 5%, 3 = 10% and so on until 21 = 100% (7 used for all experiments up to 23-05-2019)
Dp = 60.0*100.0      # Pressure layer above the surface for which the HDW index is calculated (Pa)
                     # Srock et al. use 500m, ~60 hPa at sea level.  Also grid spacing is 50 hPa above 900 hPa

if PorT == 'P':
  print 'Plume top defined by pressure minimum: ',Pmin/100.0, 'hPa'
elif PorT == 'T':
  print 'Plume top defined by Temperature minumum: ',Tmin - 273.15, 'deg C'
else:
  print '******* WARNING: No plume top temperature or pressure has been set.'
  print 'Check Pmin, Tmin and PorT'

# Dimensions
nt   = len(ncfile.dimensions['time']) # nt=24
nlat = len(ncfile.dimensions['lat'])
nlon = len(ncfile.dimensions['lon'])
nlvl = len(ncfile.dimensions['lvl'])
print(nt)

# Coordinates
t    = ncfile.variables['time'][:]
lat  = ncfile.variables['lat' ][:]
lon  = ncfile.variables['lon' ][:]
lvl  = ncfile.variables['lvl' ][:]

# Model run date and time
fd = ncfile.variables['base_date'][0]
hr = ncfile.variables['base_time'][0]
fdate = str(fd)
if hr == 0:
  hour = '0' + str(hr) + '00'
elif hr == 600:
  hour = '0' + str(hr)
else:
  hour = str(hr)

# If the time units are hours, convert to seconds
print 'Before: t[0]',t[0]
if t[0] < 20:
  t = t*24*3600.0
print 'After: t[0]',t[0]

# Find the month
V1 = round(fd/100.0)
V2 = round(V1/100.0)
V3 = V2*100.0
Month = V1 - V3
print 'Month',Month,V1, V2, V3

#st = int(round(t[0]/3600.0))
#if st < 10:
#  step = '  ' + str(st)
#elif st < 100:
#  step = ' ' + str(st)
#else:
#  step = str(st)

print fdate
print hour
  

# Determine first and last time indices
t1 = 0 
t2 = nt
#t2 = 2

print t1
print t2

for ix in range(0,nlon):
   if lon[ix] <= lonmin and lon[ix+1] >= lonmin:
     lon1 = ix
   if lon[ix] <= lonmax and lon[ix+1] >= lonmax:
     lon2 = ix+1
     break

print lon1,lon2
print lon[lon1],lon[lon2]

for ix in range(0,nlat):
   if lat[ix] >= latmax and lat[ix+1] <= latmax:
     lat2 = ix
   if lat[ix] >= latmin and lat[ix+1] <= latmin:
     lat1 = ix+1
     break

print lat1,lat2
print lat[lat1],lat[lat2]

xl = (lon[lon1], lon[lon2-1]) 
yl = (lat[lat1-1], lat[lat2])

# Load pressure level variable arrays
TT1 = ncfile.variables['air_temp'][:,:,:,:]       # Units K
RH1 = ncfile.variables['relhum'][:,:,:,:]         # Units %
uu1 = ncfile.variables['zonal_wnd'][:,:,:,:]      # Units m/s
vv1 = ncfile.variables['merid_wnd'][:,:,:,:]      # Units m/s

vars = ncfile.variables.keys()
if 'vertical_wnd' in vars:
  print 'File contains w'
  ww1 = ncfile.variables['vertical_wnd'][:,:,:,:]   # Units m/s
elif 'omega' in vars:
  print 'File contains omega'
  om1 = ncfile.variables['omega'][:,:,:,:]          # Units Pa/s
  pr4 = lvl[np.newaxis,:,np.newaxis,np.newaxis]
  ww1 = om1*0.0
  ww1 = om1*Rd*TT1/(grav*pr4*100.0)
else:
   print 'No vertical velocity variable found'


# Load surface variable arrays
zsfc1 = ncfile3.variables['topog'][:,:,:]         # Units m   Changed from ncfile2, so that zsfc1 now comes from the analysis file (09-09-2019)

psfc1 = ncfile2.variables['sfc_pres'][:,:,:]      # Units Pa
Tsfc1 = ncfile2.variables['temp_scrn'][:,:,:]     # Units K
Usfc1 = ncfile2.variables['u10'][:,:,:]          # Units m/s
Vsfc1 = ncfile2.variables['v10'][:,:,:]          # Units m/s
qsfc1 = ncfile2.variables['qsair_scrn'][:,:,:]    # Units kg/kg

# Declare additional single-column arrays
qq = lvl*0.0               # Specific Humidity. Units kg/kg
th = lvl*0.0               # Potential temperature. Units K
pr = lvl*100.0             # Pressure. Units Pa
TT = lvl*0.0               # Air Temperature. Units K
uu = lvl*0.0               # Zonal wind. Units m/s
vv = lvl*0.0               # Meridional wind. Units m/s
ww = lvl*0.0               # Vertical wind. Units m/s
dp = lvl*0.0               # Dew point temperature.  Units K

# Declare additional horizontal arrays
Hflux1 = psfc1*0.0       # Total heat flux. Units W
HfluxF = psfc1*0.0-1000. # Hflux divided by 10 m wind speed.  Hflux "flag". Based on assumption that fire intensity is proportional to wind speed.
#Bflux1 = psfc1*0.0       # Total buoyancy flux. Units m^4.s^-3
zFC1 = psfc1*0.0         # Free-convection height.  Units m
#pFC1 = psfc1*0.0         # Free-convection pressure. Units
betaFC1 = psfc1*0.0      # Free-convection buoyancy. Units m/s^2
DthFC1 = psfc1*0.0       # Free-convection plume excess potential temperature.  Units K
UML1 = psfc1*0.0         # Mean horizontal windspeed. Units m/s
Um1 = psfc1*0.0          # Mean zonal wind. Units m/s
Vm1 = psfc1*0.0          # Mean meridional wind. Units m/s
CHI = psfc1*0.0 - 5.0    # Array for C-Haines index (initialise with a value of -5, so that ocean values are off the colour scale)
FFDI= psfc1*0.0 - 15.0   # Array for Forest Fire Danger Index (assuming Drought Factor = 10)
HDW = psfc1*0.0 - 150.0  # Hot Dry Windy index (Srock et al. 2018, Atmosphere)
CHIf = psfc1*0.0 - 5.0   # Array for filtered CHI
FFDIf= psfc1*0.0 - 15.0  # Array for filtered FFDI
HDWf= psfc1*0.0 - 150.0  # Array for filtered HDW
VSTA= psfc1*0.0 - 150.0  # Array for Vesta function (tested for a PFT flag)
Pmin1=psfc1*0.0          # Array containing pressure values at the height T = Tmin.  Used if PorT set to 'T'

#  WARNING - It is important to multiply lvl/psfc1 above by some value, otherwise the
#            declared arrays will always take the same value (e.g., changing qq changes th etc.)

#  Set plotting parameters here, so that the maximum values can be used to define ocean points
# Hflux
Hf_int = 100.0
Hf_min = 0.0
Hf_max = 1100.0
Hf_con = np.arange(Hf_min,Hf_max,Hf_int)
hmap = plt.cm.get_cmap('YlOrRd_r') # Get colormap

# Hflux flag
HfF_int = 0.5
HfF_min = 0.0
HfF_max = 5.5
HfF_con = np.arange(HfF_min,HfF_max,HfF_int)
hFmap = plt.cm.get_cmap('PuBuGn_r') # Get colormap

# HDW index Contour
HDW_int = 200.0
HDW_min = 0.0
HDW_max = 790.0
HDW_con1 = np.arange(HDW_min,HDW_max,HDW_int)

HDW_int = 100.0
HDW_min = 0.0
HDW_max = 1100.0
HDW_con2 = np.arange(HDW_min,HDW_max,HDW_int)

HDW_single_con = 800.0
#HDW_con3 = np.arange(200.0,600.0,200.0)

HDWmap = plt.cm.get_cmap('PuBuGn') # Get colormap

# FFDI Contour
FFDI_int = 12.5
FFDI_min = 0.0
FFDI_max = 100.0
FFDI_con = np.arange(FFDI_min,FFDI_max,FFDI_int)
FFDImap = plt.cm.get_cmap('YlOrRd') # Get colormap

FFDI_single_con = 50.0
#FFDI_con3 = np.arange(25.0,75.0,25.0)

# Vesta Contour
VSTA_int = 2.5
VSTA_min = 0.0
VSTA_max = 32.5
VSTA_con = np.arange(VSTA_min,VSTA_max,VSTA_int)
VSTAmap = plt.cm.get_cmap('YlOrRd') # Get colormap

VSTA_single_con = 15.0


# CHI Contour
CHI_int = 2.0
CHI_min = -2.0
CHI_max = 20.0
CHI_con = np.arange(CHI_min,CHI_max,CHI_int)
CHImap = plt.cm.get_cmap('YlGnBu') # Get colormap
 
CHI_single_con = 10.0
#CHI_con3 = np.arange(4.0,10.0,2.0)

# Mixed layer windspeed
S_int = 2.0
S_min = 0.0
S_max = 22.0
scon = np.arange(S_min,S_max,S_int) # Plotting interval for wind speed
smap = plt.cm.get_cmap('YlGnBu_r') # Get colormap

# Free-convection buoyancy
BFC_int = 0.05
BFC_min = 0.0
BFC_max = 0.55
BFC_con = np.arange(BFC_min,BFC_max,BFC_int)
bmap = plt.cm.get_cmap('PuBuGn_r') # Get colormap

# Free-convection plume excess potential temeprature
DTHFC_int = 1.5
DTHFC_min = 0.0
DTHFC_max = 16.5
DTHFC_con = np.arange(DTHFC_min,DTHFC_max,DTHFC_int)
dthmap = plt.cm.get_cmap('PuBuGn_r') # Get colormap


# Free-convection height
zFC_int = 500.0
zFC_min = 0.0
zFC_max = 5500.0
zFC_con = np.arange(zFC_min,zFC_max,zFC_int)
zmap = plt.cm.get_cmap('YlOrBr_r') # Get colormap

# Determine the vector plotting frequency
nln = lon2 - lon1
nlt = lat1 - lat2
qx = int(round(max(nln,nlt)/24.0))
print 'qx',qx

#x=1.0
#print np.exp(x)
#print TT1.shape
#print th.shape
#print lvl.shape
#print lon.shape

#for iz in range(0,nlvl):
#   print iz,lvl[iz]
#print lon[450],lon[840]
#print lat[426],lat[518]
#print lon[768]
#print lat[444]

# Find pressure level at which T = Tmin
if PorT == 'T':
  for it in range(t1,t2):
    for ix in range(lon1,lon2):
      for iy in range(lat2,lat1):
       if zsfc1[0,iy,ix] != 0.0:                   # Changed from zsfc1[it,iy,ix], because zsfc1 now only has one time dimension (09-09-2019)
        for iz in range(nlvl-1,-1,-1):
#          print 'P,T',pr[iz],TT1[it,iz,iy,ix],iz
          if TT1[it,iz,iy,ix] >= Tmin:
            zdn = iz
            alpha = (Tmin - TT1[it,zdn+1,iy,ix])/(TT1[it,zdn,iy,ix] - TT1[it,zdn+1,iy,ix])
            Pmin1[it,iy,ix] = alpha*( pr[zdn] - pr[zdn+1] ) + pr[zdn+1]
#            print 'alpha,Tmin,TT1[it,zdn+1,iy,ix],TT1[it,zdn,iy,ix]',alpha,Tmin,TT1[it,zdn+1,iy,ix],TT1[it,zdn,iy,ix]
#            print 'Pmin1,pr[zdn+1],pr[zdn]',Pmin1[it,iy,ix],pr[zdn+1],pr[zdn]
            break
else:
  Pmin1 = Pmin
      

# Find 850 and 700 hPa pressure level indices for C-Haines calculation
k850 = 0
k700 = 0
for iz in range(0,nlvl):
  if pr[iz] == 85000.0:
     k850 = iz
  if pr[iz] == 70000.0:
     k700 = iz

print 'k850 = ',k850, pr[k850]
print 'k700 = ',k700, pr[k700]

# Loop through each horizontal grid point
#it=0
if (Month >= 4) and (Month <= 9):
  if (lat[lat1] <= 0.0):
    VSTS = 1 # Winter SH
  else:
    VSTS = 2 # Summer NH
else:
  if (lat[lat1] <= 0.0):
    VSTS = 2 # Summer SH
  else:
    VSTS = 1 # Winter NH

for it in range(t1,t2):
  hourUTC = np.mod(hr + t[it]/36.0,2400)
  print 'hourUTC',hourUTC

  for ix in range(lon1,lon2):
     hourOffset = 100.0*lon[ix]/15.0
     hourlocal = hourUTC+hourOffset
     if hourlocal > 2400:
       hourlocal = hourlocal - 2400
     if hourlocal < 0:
       hourlocal = hourlocal + 2400
#     print 'hourOffset,hourlocal',hourOffset,hourlocal
     if (hourlocal >= 0) and (hourlocal < 600):
       VSTP = 3
     elif (hourlocal >= 600) and (hourlocal < 1200):
       VSTP = 2
     elif (hourlocal >= 1200) and (hourlocal <= 1700):
       if VSTS == 2:
         VSTP = 1
       else:
         VSTP = 2
     elif (hourlocal >= 1700) and (hourlocal <= 1800):
       VSTP = 2
     else:
       VSTP = 3
     print 'Vesta Period',VSTP,hourlocal,VSTS,lon[ix],ix
     for iy in range(lat2,lat1):
#        print iy
#  for ix in range (749,750):
#     print ix
#     for iy in range (500,501):
#        print iy
#        print lat[iy],lon[ix]
        # Set surface variables
        zsfc = zsfc1[0,iy,ix]    # Changed from zsfc1[it,iy,ix] 09-09-2019
        psfc = psfc1[it,iy,ix]
        Tsfc = Tsfc1[it,iy,ix]
        qsfc = qsfc1[it,iy,ix]
        usfc = 0.5*( Usfc1[it,iy,ix] + Usfc1[it,iy,ix-1] )  # Convert from staggered grid to standard grid
        vsfc = 0.5*( Vsfc1[it,iy,ix] + Vsfc1[it,iy-1,ix] )
        Pmin2 = Pmin1[it,iy,ix]

#        print 'usfc,vsfc ',usfc,vsfc,ix,iy

#        # Initialise HDW maxima variables
#        VPDmax = 0.0
#        SPDmax = 0.0

        # Calculate specific humidity and potential temperature
        for iz in range(0,nlvl):
           #  saturation vapour pressure
           e_star = 611.2*np.exp( (17.67*TT1[it,iz,iy,ix]- 4826.56)/(TT1[it,iz,iy,ix] - 29.65) )
           #  vapour pressure
           ee = RH1[it,iz,iy,ix]*e_star/100.0
           #  specific humidity
           qq[iz] = (epsln*ee)/( pr[iz] - ee*(1.0-epsln) )
           #  mixing ratio
           rr = qq[iz]/(1.0 - qq[iz])
           #  potential temperature
           th[iz] = TT1[it,iz,iy,ix]*(p0/pr[iz])**(Rd/Cpd*(1 - 0.24*rr))
           # Dew point temperature
           dp[iz] = 243.5/( ( 17.67/np.log(ee/611.2) ) - 1.0) + 273.15
           # Fill remaining single column arrays
           TT[iz] = TT1[it,iz,iy,ix]
           uu[iz] = uu1[it,iz,iy,ix]
           vv[iz] = vv1[it,iz,iy,ix]
           ww[iz] = ww1[it,iz,iy,ix]
#           print 'ww[iz] ',ww[iz],iz
#           print pr[iz]*0.01,TT[iz],TT[iz]-273.15,qq[iz]*1000.0,iz

#           # Begin HDW calculation, but only over land (zsfc > 0)
#           if zsfc != 0:
#             if pr[iz] < psfc and pr[iz] > psfc-Dp:
#               Ts = th[iz]*(psfc/p0)**(Rd/Cpd*(1 - 0.24*rr))
#               e_star_sfc = 611.2*np.exp( (17.67*Ts- 4826.56)/(Ts - 29.65) )
#               e_sfc = qq[iz]*psfc/( epsln + qq[iz]*(1.0 - epsln) )
#               VPD = e_star_sfc - e_sfc
##               print psfc,Dp,rr,qq[iz],pr[iz]
##               print Ts,e_star_sfc,e_sfc,VPD
#               if VPD > VPDmax:
#                  VPDmax = VPD
#               SPD = (uu[iz]**2 + vv[iz]**2)**0.5
##               print SPD, uu[iz], vv[iz]
#               if SPD > SPDmax:
#                  SPDmax = SPD
#             
#           ### - End of HDW if construct
#        ### - End of iz loop
##        print SPDmax, VPDmax

        # Call Fortran module to calculate fire-power etc. 
#        print 'before heat_flux_calc'
#        print phi
#        print DbSP
#        print ni,nj
#        print zp
#        print beta_e
#        print Pmin
#        print betaSPmin,Wmin,Umin,Prcntg
        if zsfc != 0:
          [UML,Um,Vm,betaFC,DthFC,zFC,pFC,Bflux,Hflux]=subroutines.heat_flux_calc(TT,qq,uu,vv,ww,th,pr,nlvl,zsfc,psfc,Tsfc,\
                                                     phi,DbSP,ni,nj,zp,beta_e,Pmin2,betaSPmin,Wmin,Umin,Prcntg)
        else:
          UML=S_max*1.0
          Um =0.0
          Vm =0.0
          betaFC=BFC_max*1.0
          DthFC=DTHFC_max*1.0
          zFC=zFC_max*1.0
          pFC=100.0
          Bflux=1.0e9
          Hflux=2*Hf_max*1.0e9
         # End of if construct

#        print UML,Um,Vm,betaFC,zFC,pFC,Bflux,Hflux

        # Transfer heat_flux_calc output to arrays
        UML1[it,iy,ix] = UML
        Um1[it,iy,ix] = Um
        Vm1[it,iy,ix] = Vm
        betaFC1[it,iy,ix] = betaFC
        DthFC1[it,iy,ix] = DthFC
        zFC1[it,iy,ix] = zFC
#        pFC1[it,iy,ix] = pFC
#        Bflux1[it,iy,ix] = Bflux
        Hflux1[it,iy,ix] = Hflux    

#        # Fill HDW array 
#        HDW[it,iy,ix] = VPDmax*SPDmax

        # Calulate C-Haines index, FFDI, Vesta and Hflux Flagged for land points
        if zsfc != 0:
          DD850 = min(TT[k850] - dp[k850],30.0)
          CA = 0.5*(TT[k850] - TT[k700] - 4)
          CB = DD850*0.3333 - 1.0 
          CHI[it,iy,ix] = CA + CB
#          print 'C-Haines check'
#          print 'TT850, DP850, TT700 ',TT[k850],dp[k850],TT[k700]
#          print 'CA, CB, CHI ',CA, CB, CHI[it,iy,ix]
          e_star_sfc = 611.2*np.exp( (17.67*Tsfc- 4826.56)/(Tsfc - 29.65) )
          e_sfc = qsfc*psfc/( epsln + qsfc*(1.0 - epsln) )
          RHsfc = (e_sfc/e_star_sfc)*100.0
          spdsfc = (usfc**2 + vsfc**2)**0.5
#          FFDI[it,iy,ix] = 2.0*np.exp(1.82265 - 0.0345*RHsfc + 0.0338*(Tsfc-273.15) + 0.08424*spdsfc ) # Drought Factor = 10
          VSTA_V = 1.5*(max(spdsfc,3.0)-1.39)**0.858
          if VSTP == 1:
            VSTA_M = 2.76 + 0.124*RHsfc - 0.0187*(Tsfc-273.15)  # Period 1
          elif VSTP == 2:
            VSTA_M = 3.60 + 0.169*RHsfc - 0.0450*(Tsfc-273.15)  # Period 2
          elif VSTP == 3:
            VSTA_M = 3.08 + 0.198*RHsfc - 0.0483*(Tsfc-273.15)  # Period 3
#          VSTA_M = 2.76 + 0.124*RHsfc - 0.0187*(Tsfc-273.15)  # Period 1 To be used alone for arctic mid-summer
          VSTA[it,iy,ix] = 18.35*(VSTA_M**(-1.495))*VSTA_V
#          print 'VSTA_V,VSTA_M,VSTA[it,iy,ix]',VSTA_V,VSTA_M,VSTA[it,iy,ix]
#          HfluxF[it,iy,ix] = Hflux/VSTA[it,iy,ix]
#          HfluxF[it,iy,ix] = Hflux/(VSTA[it,iy,ix]**2)
          if (VSTA[it,iy,ix] > 2.0) and (spdsfc > 2.0):
#            HfluxF[it,iy,ix] = Hflux/(VSTA[it,iy,ix]**2)
            HfluxF[it,iy,ix] = Hflux/(VSTA[it,iy,ix]*spdsfc)
          else:
            HfluxF[it,iy,ix] = 1.0e12

#          spdsfc = max(0.1,spdsfc)
##          HfluxF[it,iy,ix] = Hflux/spdsfc
##          HfluxF[it,iy,ix] = Hflux/(spdsfc**2)
#          if FFDI[it,iy,ix] <= 12.0:
#            HfluxF[it,iy,ix] = 1.0e12
#          else:
##            HfluxF[it,iy,ix] = Hflux/FFDI[it,iy,ix]
##            HfluxF[it,iy,ix] = Hflux/(FFDI[it,iy,ix]**2)
#            HfluxF[it,iy,ix] = Hflux/(FFDI[it,iy,ix]*spdsfc)

#          print 'FFDI check'
#          print 'RH, T, Spd ',RHsfc,Tsfc-273.15,spdsfc
#          print 'FFDI ',FFDI[it,iy,ix]
#          print 'HfluxF check', HfluxF[it,iy,ix]


###   - End of ix,iy loops 
### - End of it loop

ncfile.close()
ncfile2.close()

# Set Hflux units to giga Watts
Hflux1=Hflux1*1.0e-9
HfluxF=HfluxF*1.0e-9
# Convert from beta_fc to b_fc
betaFC1=betaFC1*grav
## Convert HDW from units of m/s.Pa to m/s.hPa
#HDW=HDW*0.01

### - Begin plotting

### - Output file 1

for it in range(t1,t2):
   st = int(round(t[it]/3600.0))
   if st < 10:
     step = '00' + str(st)
   elif st < 100:
     step = '0' + str(st)
   else:
     step = str(st)
   print step
   hourUTC = np.mod(hr + t[it]/36.0,2400)
   ihourUTC = int(hourUTC)
   if ihourUTC < 10:
     strhourUTC = '000' + str(ihourUTC)
   elif ihourUTC < 100:
     strhourUTC = '00' + str(ihourUTC)
   elif ihourUTC < 1000:
     strhourUTC = '0' + str(ihourUTC)
   else:
     strhourUTC = str(ihourUTC)
   
   print 'strhourUTC',strhourUTC
 

   # Set spaces between sub plots
#   plt.subplots_adjust(hspace = 0.1, wspace = 0.1)

   # Create figure
   # Heat flux and winds
   fig = plt.figure(1, figsize=(16,14))
   plt.subplot(2,2,1,aspect='equal')
#   plt.contourf(lon,lat,Hflux1[it,:,:],Hf_con,cmap=hmap) # Draw filled contours
#   plt.contourf(lon,lat,Hflux1[it,:,:],levels=[00.0,25.0,50.0,100.0,150.0,225.0,300.0,375.0,450.0,525.0,600.0,700.0,800.0],cmap=hmap)
   clev = 2**np.arange(4,10.5,0.5)
   plt.contourf(lon,lat,Hflux1[it,:,:],clev,norm=col.LogNorm(),cmap=hmap)
#   cbar = plt.colorbar(format=tick.LogFormatterSciNotation(base=2.0)) # Add colorbar to the plot
   cbar = plt.colorbar(format=tick.LogFormatter(base=2.0)) # Add colorbar to the plot
#   cbar = plt.colorbar()
   # Access to cbar tick labels (change font, etc.)
   cbar.ax.tick_params(labelsize=8)
#   plt.contour(lon,lat,HDWf[it,:,:],HDW_con1,colors=['cyan','fuchsia','red'])
#   plt.contourf(lon,lat,HDW[it,:,:],HDW_con2,colors='none',hatches=['/'])
   # Plot wind vectors
   Q = plt.quiver(lon[::qx],lat[::qx],Um1[it,::qx,::qx],Vm1[it,::qx,::qx],scale=0.5/0.01,units='inches',color='k') # units='y')
   plt.quiverkey(Q, 0.1, 1.05, 10, r'$10 \frac{m}{s}$', labelpos='W', fontproperties={'weight': 'bold'})
   # Plot map of Australia
   plt.plot(map_lon,map_lat,'k-',linewidth=2.0)
   plt.xlabel('longitude', fontsize=14) # x-axis label
   plt.ylabel('latitude', fontsize=14)  # y-axis label
   plt.xlim(xl); plt.ylim(yl) # Axis limits
# Use for analysis
#   plt.title('Heat flux (GW)  ' + fdate + ': ' +hour,fontsize=10)
# Use for forecasts
   plt.title('Heat flux (GW)  ' + fdate + ': ' +hour+ ' + ' +step+ '  ['+strhourUTC+' UTC]',fontsize=14)

# Second plot
   plt.subplot(2,2,2,aspect='equal') # One subplot
   plt.contourf(lon,lat,zFC1[it,:,:],zFC_con,cmap=zmap) # draw filled contours for free-convection height
   cbar = plt.colorbar() # add colorbar to the plot
   cbar.ax.tick_params(labelsize=8)
   Q = plt.quiver(lon[::qx],lat[::qx],Um1[it,::qx,::qx],Vm1[it,::qx,::qx],scale=0.5/0.01,units='inches',color='k') # units='y')
   plt.quiverkey(Q, 0.1, 1.05, 10, r'$10 \frac{m}{s}$', labelpos='W', fontproperties={'weight': 'bold'})
   plt.plot(map_lon,map_lat,'k-',linewidth=2.0)
   plt.xlabel('longitude', fontsize=14) # x-axis label
   plt.ylabel('latitude', fontsize=14)  # y-axis label
   plt.xlim(xl); plt.ylim(yl) # axis limits
   plt.title('Free-convection height (m)  ',fontsize=14)

# Third plot
   plt.subplot(2,2,3,aspect='equal')
   plt.contourf(lon,lat,UML1[it,:,:],scon,cmap=smap) # Draw filled contours
   cbar = plt.colorbar() # Add colorbar to the plot
   # Access to cbar tick labels (change font, etc.)
   cbar.ax.tick_params(labelsize=8)
   # Plot wind vectors
   Q = plt.quiver(lon[::qx],lat[::qx],Um1[it,::qx,::qx],Vm1[it,::qx,::qx],scale=0.5/0.01,units='inches',color='k') # units='y')
   plt.quiverkey(Q, 0.1, 1.05, 10, r'$10 \frac{m}{s}$', labelpos='W', fontproperties={'weight': 'bold'})
   # Plot map of Australia
   plt.plot(map_lon,map_lat,'k-',linewidth=2.0)
   plt.xlabel('longitude', fontsize=14) # x-axis label
   plt.ylabel('latitude', fontsize=14)  # y-axis label
   plt.xlim(xl); plt.ylim(yl) # Axis limits
   plt.title('Mixed-layer wind speed (m/s)  ',fontsize=14)

# Fourth plot
   plt.subplot(2,2,4,aspect='equal') # One subplot
#   plt.contourf(lon,lat,betaFC1[it,:,:],BFC_con,cmap=bmap) # draw filled contours
   plt.contourf(lon,lat,DthFC1[it,:,:],DTHFC_con,cmap=dthmap) # draw filled contours
   cbar = plt.colorbar() # add colorbar to the plot
   cbar.ax.tick_params(labelsize=8)
   Q = plt.quiver(lon[::qx],lat[::qx],Um1[it,::qx,::qx],Vm1[it,::qx,::qx],scale=0.5/0.01,units='inches',color='k') # units='y')
   plt.quiverkey(Q, 0.1, 1.05, 10, r'$10 \frac{m}{s}$', labelpos='W', fontproperties={'weight': 'bold'})
   plt.plot(map_lon,map_lat,'k-',linewidth=2.0)
   plt.xlabel('longitude', fontsize=14) # x-axis label
   plt.ylabel('latitude', fontsize=14)  # y-axis label
   plt.xlim(xl); plt.ylim(yl) # axis limits
   plt.title('Free-convection Delta theta (K)  ',fontsize=14)

   plt.tight_layout(h_pad=1)
   # Add warning label in the middle of the figure
   plt.suptitle('WARNING: Product not operationally supported. Provided for research testing and feedback purposes only.',fontsize=20,y=0.515)

#   plt.show()
   plt.savefig(figdir + '/PFT_' + fdate + '_' +hour+ '_' +step+ '.png', dpi=100)
   plt.close(fig)

## End of it loop First file

### - Output file 2

for it in range(t1,t2):
   st = int(round(t[it]/3600.0))
   if st < 10:
     step = '00' + str(st)
   elif st < 100:
     step = '0' + str(st)
   else:
     step = str(st)
   print step
   hourUTC = np.mod(hr + t[it]/36.0,2400)
   ihourUTC = int(hourUTC)
   if ihourUTC < 10:
     strhourUTC = '000' + str(ihourUTC)
   elif ihourUTC < 100:
     strhourUTC = '00' + str(ihourUTC)
   elif ihourUTC < 1000:
     strhourUTC = '0' + str(ihourUTC)
   else:
     strhourUTC = str(ihourUTC)

   print 'strhourUTC',strhourUTC

   # Set spaces between sub plots
#   plt.subplots_adjust(hspace = 0.1, wspace = 0.1)

   # Create figure
   # Heat flux and winds
   fig = plt.figure(1, figsize=(16,14))
   plt.subplot(2,2,1,aspect='equal')
#   plt.contourf(lon,lat,Hflux1[it,:,:],Hf_con,cmap=hmap) # Draw filled contours
#   plt.contourf(lon,lat,Hflux1[it,:,:],levels=[00.0,25.0,50.0,100.0,150.0,225.0,300.0,375.0,450.0,525.0,600.0,700.0,800.0],cmap=hmap)
#   cbar = plt.colorbar() # Add colorbar to the plot
#   clev = 2**np.arange(4,10.5,0.5)
#   plt.contourf(lon,lat,Hflux1[it,:,:],clev,norm=col.LogNorm(),cmap=hmap)
#   cbar = plt.colorbar(format=tick.LogFormatter(base=2.0)) # Add colorbar to the plot
   clev = 10**np.arange(0.0,3.25,0.25)
   plt.contourf(lon,lat,Hflux1[it,:,:],clev,norm=col.LogNorm(),cmap=hmap)
   cbar = plt.colorbar(format=tick.LogFormatter(base=10.0)) # Add colorbar to the plot
   plt.contour(lon,lat,CHI[it,:,:],CHI_single_con,colors=['cyan'])
   plt.contour(lon,lat,VSTA[it,:,:],VSTA_single_con,colors=['fuchsia'])
#   plt.contour(lon,lat,HDW[it,:,:],HDW_single_con,colors=['lime'])
   # Access to cbar tick labels (change font, etc.)
   cbar.ax.tick_params(labelsize=8)
#   plt.contour(lon,lat,HDW[it,:,:],HDW_con1,colors=['cyan','fuchsia','red'])
   # Plot wind vectors
   Q = plt.quiver(lon[::qx],lat[::qx],Um1[it,::qx,::qx],Vm1[it,::qx,::qx],scale=0.5/0.01,units='inches',color='k') # units='y')
   plt.quiverkey(Q, 0.1, 1.05, 10, r'$10 \frac{m}{s}$', labelpos='W', fontproperties={'weight': 'bold'})
   # Plot map of Australia
   plt.plot(map_lon,map_lat,'k-',linewidth=2.0)
   plt.xlabel('longitude', fontsize=14) # x-axis label
   plt.ylabel('latitude', fontsize=14)  # y-axis label
   plt.xlim(xl); plt.ylim(yl) # Axis limits
# Use for analysis
#   plt.title('Heat flux (GW)  ' + fdate + ': ' +hour,fontsize=10)
# Use for forecasts
   plt.title('Heat flux (GW)  ' + fdate + ': ' +hour+ ' + ' +step+'  ['+strhourUTC+' UTC]',fontsize=14)
   # plt.title(str(lvl[lev]) + ' hPa winds, geop. height (blue) and speed (colour, m/s) ' + date + ' at ' + str(base_time[it]),fontsize=10)  # Plot title

# Second plot
   plt.subplot(2,2,2,aspect='equal') # One subplot
#   plt.contourf(lon,lat,HDW[it,:,:],HDW_con2,cmap=HDWmap) # draw filled contours
   plt.contourf(lon,lat,HfluxF[it,:,:],HfF_con,cmap=hFmap) # draw filled contours
   cbar = plt.colorbar() # add colorbar to the plot
   cbar.ax.tick_params(labelsize=8)
#   BFCcon=plt.contour(lon,lat,betaFC1[it,:,:],BFC_con,colors='b')
#   plt.clabel(BFCcon,BFCcon.levels[::2],fontsize=10)
   Q = plt.quiver(lon[::qx],lat[::qx],Um1[it,::qx,::qx],Vm1[it,::qx,::qx],scale=0.5/0.01,units='inches',color='k') # units='y')
   plt.quiverkey(Q, 0.1, 1.05, 10, r'$10 \frac{m}{s}$', labelpos='W', fontproperties={'weight': 'bold'})
   plt.plot(map_lon,map_lat,'k-',linewidth=2.0)
   plt.xlabel('longitude', fontsize=14) # x-axis label
   plt.ylabel('latitude', fontsize=14)  # y-axis label
   plt.xlim(xl); plt.ylim(yl) # axis limits
#   plt.title('Hot Dry Windy Index (m/s.hPa)  ',fontsize=14)
#   plt.title('PFT Flag ( GW/(m/s)^2 )  ',fontsize=14)
#   plt.title('PFT Flag ( PFT/FFDI )  ',fontsize=14)
#   plt.title('PFT Flag ( PFT/(FFDI)^2 )  ',fontsize=14)
#   plt.title('PFT Flag ( PFT/(FFDI*Uspd) )  ',fontsize=14)
#   plt.title('PFT Flag ( PFT/Vesta ) )  ',fontsize=14)
#   plt.title('PFT Flag ( PFT/Vesta**2 ) )  ',fontsize=14)
   plt.title('PFT Flag ( PFT/(Vesta*Uspd) )  ',fontsize=14)

# Third plot
   plt.subplot(2,2,3,aspect='equal')
   plt.contourf(lon,lat,VSTA[it,:,:],VSTA_con,cmap=VSTAmap) # Draw filled contours
#   plt.contourf(lon,lat,FFDI[it,:,:],FFDI_con,cmap=FFDImap) # Draw filled contours
#   plt.contourf(lon,lat,FFDI[it,:,:],colors=['darkseagreen','skyblue','yellow','orange','r','firebrick'],cmap=FFDImap)
   cbar = plt.colorbar() # Add colorbar to the plot
   plt.contour(lon,lat,VSTA[it,:,:],VSTA_single_con,colors=['fuchsia'])
   # Access to cbar tick labels (change font, etc.)
   cbar.ax.tick_params(labelsize=8)
   # Plot wind vectors
   Q = plt.quiver(lon[::qx],lat[::qx],Um1[it,::qx,::qx],Vm1[it,::qx,::qx],scale=0.5/0.01,units='inches',color='k') # units='y')
   plt.quiverkey(Q, 0.1, 1.05, 10, r'$10 \frac{m}{s}$', labelpos='W', fontproperties={'weight': 'bold'})
   # Plot map of Australia
   plt.plot(map_lon,map_lat,'k-',linewidth=2.0)
   plt.xlabel('longitude', fontsize=14) # x-axis label
   plt.ylabel('latitude', fontsize=14)  # y-axis label
   plt.xlim(xl); plt.ylim(yl) # Axis limits
#   plt.title('Forest Fire Danger Index  ',fontsize=14)
   plt.title('Vesta function ',fontsize=14)

# Fourth plot
   plt.subplot(2,2,4,aspect='equal') # One subplot
   plt.contourf(lon,lat,CHI[it,:,:],CHI_con,cmap=CHImap) # draw filled contours for C-Haines index
   cbar = plt.colorbar() # add colorbar to the plot
   plt.contour(lon,lat,CHI[it,:,:],CHI_single_con,colors=['cyan'])
   cbar.ax.tick_params(labelsize=8)
   Q = plt.quiver(lon[::qx],lat[::qx],Um1[it,::qx,::qx],Vm1[it,::qx,::qx],scale=0.5/0.01,units='inches',color='k') # units='y')
   plt.quiverkey(Q, 0.1, 1.05, 10, r'$10 \frac{m}{s}$', labelpos='W', fontproperties={'weight': 'bold'})
   plt.plot(map_lon,map_lat,'k-',linewidth=2.0)
   plt.xlabel('longitude', fontsize=14) # x-axis label
   plt.ylabel('latitude', fontsize=14)  # y-axis label
   plt.xlim(xl); plt.ylim(yl) # axis limits
   plt.title('C-Haines ',fontsize=14)

   plt.tight_layout(h_pad=1)

   # Add warning label in the middle of the figure
   plt.suptitle('WARNING: Product not operationally supported. Provided for research testing and feedback purposes only.',fontsize=20,y=0.515)

#   plt.show()
   plt.savefig(figdir + '/PFT2_' + fdate + '_' +hour+ '_' +step+ '.png', dpi=100)
   plt.close(fig)

## End of it loop Second file


print 'End of script'
