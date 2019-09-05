# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 09:05:52 2016

@author: Jeff
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter#, MultipleLocator
from matplotlib.projections import register_projection
import numpy as np

import skewt
import thermo

def plot_dryadiabat(ax,theta,p,p0):
    # p should be an array of pressures at which the adiabat is plotted
    T = theta*(p/p0)**thermo.kappa
    ax.semilogy(T, p, 'c-')
    
def plot_mixrat(ax,r,p):
    # p should be an array of pressures (hPa) at which the adiabat is plotted
    e = p * r / (thermo.eps + r)
    loge = np.log(e)
    Td = (243.5*loge - 440.8)/(19.48 - loge) + 273.15
    ax.semilogy(Td,p,'c-')
 
    
mpl.rcParams['font.size'] = 18.0
mpl.rcParams['mathtext.fontset'] = 'stixsans'
mtfs = 22.0  # fontsize for math text

register_projection(skewt.SkewXAxes)

p0 = 1e5

th_bl = 303
q_bl = 0.005

T_bl = th_bl*(p0/1e5)**(thermo.kappa)
r_bl = q_bl/(1 - q_bl)
e_bl = p0 * r_bl / (thermo.eps + r_bl)
loge_bl = np.log(0.01*e_bl)
Td_bl = (243.5*loge_bl - 440.8)/(19.48 - loge_bl) + 273.15

gamma = 2.0
delta = 6.6*1e3  # need the 1e3 to get units in K/(kg/kg)

th_fire = gamma*th_bl
T_fire = th_fire*(p0/1e5)**(thermo.kappa)
q_fire = ((gamma-1)/delta)*th_bl + 0.86*q_bl
r_fire = q_fire/(1 - q_fire)
e_fire = p0 * r_fire / (thermo.eps + r_fire)
loge_fire = np.log(0.01*e_fire)
Td_fire = (243.5*loge_fire - 440.8)/(19.48 - loge_fire) + 273.15

al = np.linspace(0.0,1.0,1000)

th_sp = (1-al)*th_bl + al*th_fire
q_sp = (1-al)*q_bl + al*q_fire
r_sp = q_sp/(1 - q_sp)

T0_sp = th_sp*(p0/1e5)**(thermo.Rd/thermo.Cpd)

e = p0 * r_sp / (thermo.eps + r_sp)
loge = np.log( 0.01*np.maximum(e, 1e-100) )
T_sp = 2840.0 / (3.5*np.log(T0_sp) - loge - 4.805) + 55.0
p_sp = 1e5*(T_sp/th_sp)**(thermo.Cpd/thermo.Rd)

fig = plt.figure(1)#figsize=(6.5875, 6.2125))
ax = fig.add_subplot(111, projection='skewx')

yticks = np.ravel(np.outer(10**np.arange(0,3),np.linspace(1,10,10)))

plt.grid(True)
for thx in range(300,2000,100):
    plot_dryadiabat(ax,thx,yticks,1e3)
for rx in [1e-4,2e-4,4e-4,1e-3,2e-3,4e-3,1e-2,2e-2,4e-2,0.1,0.2,0.4,1.0]:
    plot_mixrat(ax,rx,yticks)

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dicatated by the typical meteorological plot
ax.semilogy(T_sp, p_sp*1e-2, 'r')
#ax.semilogy([T_fire,T_sp[-1],Td_fire],[p0*1e-2,p_sp[-1]*1e-2,p0*1e-2],'b-')
ax.semilogy([T_bl,T_sp[0],Td_bl],[p0*1e-2,p_sp[0]*1e-2,p0*1e-2],'b-')

# Disables the log-formatting that comes with semilogy
ax.yaxis.set_major_formatter(ScalarFormatter())
ax.set_yticks(yticks)
ax.set_ylim(1050, 50)

#ax.xaxis.set_major_locator(MultipleLocator(10))
ax.set_xlim(Td_bl,400)

ax.text(0.8,0.9,  r'$\theta_{bl}='  +'{:.0f}'.format(th_bl)   +'\mathrm{K}$',transform=ax.transAxes,fontsize=mtfs)
ax.text(0.8,0.825,r'$q_{bl}='       +'{:.1f}'.format(q_bl*1e3)+'\mathrm{g/kg}$',transform=ax.transAxes,fontsize=mtfs)
ax.text(0.8,0.75, r'$\gamma='       +'{:.1f}$'.format(gamma),transform=ax.transAxes,fontsize=mtfs)
ax.text(0.8,0.675,r'$\delta='       +'{:.1f}$'.format(delta*1e-3),transform=ax.transAxes,fontsize=mtfs)
ax.text(0.8,0.6,  r'$\theta_{fire}='+'{:.0f}'.format(th_fire)+'\mathrm{K}$',transform=ax.transAxes,fontsize=mtfs)
ax.text(0.8,0.525,r'$q_{fire}='     +'{:.1f}'.format(q_fire*1e3)+'\mathrm{g/kg}$',transform=ax.transAxes,fontsize=mtfs)
ax.set_xlabel('$T \, (\mathrm{K})$')
ax.set_ylabel('$p \, (\mathrm{hPa})$')

fig.set_size_inches(20,10)
#plt.show()
plt.savefig('skewT_th{:.0f}_q{:.0f}_gam{:.0f}_del{:.0f}.png'.format(th_bl,q_bl*1e3,gamma,delta*1e-3),dpi=200)
    
    