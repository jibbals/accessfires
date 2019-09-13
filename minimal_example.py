# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 15:50:42 2019
    minimal example skewT problem
@author: Jesse
"""

import matplotlib.pyplot as plt
from metpy.units import units
from metpy.plots import SkewT
import numpy as np

# Pressure in hPa
p = np.array([1006.0437475 ,  995.76678167,  971.08324364,  932.77561992,
              882.03055039,  820.47928304,  750.09469793,  673.30749624,
              593.14227853,  512.06441041,  431.84615956,  343.13325525,
              206.98255942]) * units.hPa
# Temperature in Kelvin
T = np.array([308.75318763, 306.06809246, 303.6426414 , 300.10817166,
              295.43166071, 289.50015275, 283.0458144 , 276.37844266,
              272.42841696, 264.21209355, 254.43607077, 240.42762839,
              215.55300601]) * units.degK
# Dewpoint temperature in Kelvin
Td = np.array([287.13853628, 286.57071724, 286.14560658, 285.52678341,
               284.70488872, 283.57435228, 281.99306255, 274.24633018,
               254.00050803, 231.5517735 , 225.75652085, 220.79078678,
               206.97078861]) * units.degK

# figure with kelvin on x-axis
skew = SkewT()
skew.plot(p,T,'k')
skew.plot(p,Td,'g')
skew.ax.set_xlim(240,330)
skew.plot_dry_adiabats() # 
skew.plot_moist_adiabats()
skew.plot_mixing_lines()
plt.title("F160 plot using Kelvin (wrong adiabats/mixing ratio?)")
plt.savefig('mintest1.png')
print('info: mintest figure saved')
plt.close()

# figure with celcius on x-axis
skew = SkewT()
skew.plot(p,T.to(units.degC),'k')
skew.plot(p,Td.to(units.degC),'g')
skew.ax.set_xlim(-30,60)
skew.plot_dry_adiabats() # 
skew.plot_moist_adiabats()
skew.plot_mixing_lines()
plt.title("F160 plot using Celcius")
plt.savefig('mintest2.png')
print('info: mintest figure saved')
plt.close()

skew = SkewT()
skew.plot(p,T,'k')
skew.plot(p,Td,'g')
skew.ax.set_xlim(240,330)
skew.plot_dry_adiabats(t0=np.arange(250,380,10)) # 
skew.plot_moist_adiabats(t0=np.arange(250,380,10))
skew.plot_mixing_lines()
plt.title("F160 plot using Kelvin and no units in t0")
plt.savefig('mintest3.png')
print('info: mintest figure saved')
plt.close()

skew = SkewT()
skew.plot(p,T.to(units.degC),'k')
skew.plot(p,Td.to(units.degC),'g')
skew.ax.set_xlim(-30,60)
skew.plot_dry_adiabats(t0=np.arange(-30,60,10)) # 
skew.plot_moist_adiabats(t0=np.arange(-30,61,10))
skew.plot_mixing_lines()
plt.title("F160 plot using Celcius and no units in t0")
plt.savefig('mintest4.png')
print('info: mintest figure saved')
plt.close()

# Temperature in Kelvin without units
T = np.array([308.75318763, 306.06809246, 303.6426414 , 300.10817166,
              295.43166071, 289.50015275, 283.0458144 , 276.37844266,
              272.42841696, 264.21209355, 254.43607077, 240.42762839,
              215.55300601])
# Dewpoint temperature in Kelvin without units
Td = np.array([287.13853628, 286.57071724, 286.14560658, 285.52678341,
               284.70488872, 283.57435228, 281.99306255, 274.24633018,
               254.00050803, 231.5517735 , 225.75652085, 220.79078678,
               206.97078861])

skew = SkewT()

skew.plot(p,T,'k')
print(skew.ax.get_xlim())
skew.plot(p,Td,'g')
skew.ax.set_xlim(240,330)
print(skew.ax.xaxis.units) # units are none
skew.plot_dry_adiabats() # 
print(skew.ax.get_xlim())
skew.plot_moist_adiabats()
skew.plot_mixing_lines()
plt.title("F160 plot using non-dimensional inputs")
plt.savefig('mintest5.png')
print('info: mintest figure saved')
plt.close()

#from datetime import datetime 
#from utilities import plotting, utils, fio, constants
## Use datetime and latlon to determine what data to read
#extentname='waroona'
#extent = plotting._extents_[extentname]
#    
## read pressure and temperature cubes
#_,ro1,th1,_= fio.read_waroona(datetime(2016,1,6,5), extent=extent, add_dewpoint=True, add_winds=True)#, add_theta=True)
#p,t,td = th1.extract(['air_pressure','air_temperature','dewpoint_temperature'])
#pro, u, v  = ro1.extract(['air_pressure','u','v'])
#    
#prof=(0,slice(None,130,10),50,50)
#p[prof].data.data
#t[prof].data.data
#td[prof].data.data


###################
###################
##
#def plot_dry_adiabats(self, t0=None, p=None, **kwargs):
#        r"""Plot dry adiabats.
#
#        Adds dry adiabats (lines of constant potential temperature) to the
#        plot. The default style of these lines is dashed red lines with an alpha
#        value of 0.5. These can be overridden using keyword arguments.
#
#        Parameters
#        ----------
#        t0 : array_like, optional
#            Starting temperature values in Kelvin or Celcius. 
#            If none are given, they will be generated using the current 
#            temperature range at the bottom of the plot.
#        p : array_like, optional
#            Pressure values to be included in the dry adiabats. If not
#            specified, they will be linearly distributed across the current
#            plotted pressure range.
#        kwargs
#            Other keyword arguments to pass to :class:`matplotlib.collections.LineCollection`
#
#        Returns
#        -------
#        matplotlib.collections.LineCollection
#            instance created
#
#        See Also
#        --------
#        :func:`~metpy.calc.thermo.dry_lapse`
#        :meth:`plot_moist_adiabats`
#        :class:`matplotlib.collections.LineCollection`
#
#        """
#        # Determine set of starting temps if necessary
#        if t0 is None:
#            xmin, xmax = self.ax.get_xlim()
#            t0 = np.arange(xmin, xmax + 1, 10)
#            if self.ax.xaxis.units is not None:
#                t0 = t0 * self.ax.xaxis.units
#
#        # Get pressure levels based on ylims if necessary
#        if p is None:
#            p = np.linspace(*self.ax.get_ylim()) * units.mbar # assume mbar
#        
#        # Make sure t0 has units
#        if not hasattr(t0,'units'):
#            # Assume t0 is in celcius, unless the minimum is > 150
#            tunits = units.degC
#            if np.min(t0) > 150:
#                tunits = units.degK
#            t0=t0 * tunits
#        
#        # Assemble into data for plotting
#        t = dry_lapse(p, t0[:, np.newaxis], 1000. * units.mbar).to(t0.units)
#        linedata = [np.vstack((ti, p)).T for ti in t]
#
#        # Add to plot
#        kwargs.setdefault('colors', 'r')
#        kwargs.setdefault('linestyles', 'dashed')
#        kwargs.setdefault('alpha', 0.5)
#        return self.ax.add_collection(LineCollection(linedata, **kwargs))
#
#    def plot_moist_adiabats(self, t0=None, p=None, **kwargs):
#        r"""Plot moist adiabats.
#
#        Adds saturated pseudo-adiabats (lines of constant equivalent potential
#        temperature) to the plot. The default style of these lines is dashed
#        blue lines with an alpha value of 0.5. These can be overridden using
#        keyword arguments.
#
#        Parameters
#        ----------
#        t0 : array_like, optional
#            Starting temperature values. If none are given, they will be
#            generated using the current temperature range at the bottom of
#            the plot.
#        p : array_like, optional
#            Pressure values to be included in the moist adiabats. If not
#            specified, they will be linearly distributed across the current
#            plotted pressure range.
#        kwargs
#            Other keyword arguments to pass to :class:`matplotlib.collections.LineCollection`
#
#        Returns
#        -------
#        matplotlib.collections.LineCollection
#            instance created
#
#        See Also
#        --------
#        :func:`~metpy.calc.thermo.moist_lapse`
#        :meth:`plot_dry_adiabats`
#        :class:`matplotlib.collections.LineCollection`
#
#        """
#        # Determine set of starting temps if necessary
#        if t0 is None:
#            xmin, xmax = self.ax.get_xlim()
#            t0 = np.arange(xmin, xmax + 1, 5)
#            if self.ax.xaxis.units is not None:
#                t0 = t0 * self.ax.xaxis.units
#        
#        # Get pressure levels based on ylims if necessary
#        if p is None:
#            p = np.linspace(*self.ax.get_ylim()) * units.mbar
#        
#        # Make sure t0 has units
#        if not hasattr(t0, 'units'):
#            # Assume t0 is in celcius, unless the minimum is > 150
#            tunits = units.degC
#            if np.min(t0) > 150:
#                tunits = units.degK
#            t0 = t0 * tunits
#        
#        # Assemble into data for plotting
#        t = moist_lapse(p, t0[:, np.newaxis], 1000. * units.mbar).to(t0.units)
#        linedata = [np.vstack((ti, p)).T for ti in t]
#        
#        # Add to plot
#        kwargs.setdefault('colors', 'b')
#        kwargs.setdefault('linestyles', 'dashed')
#        kwargs.setdefault('alpha', 0.5)
#        return self.ax.add_collection(LineCollection(linedata, **kwargs))
#
#    def plot_mixing_lines(self, w=None, p=None, **kwargs):
#        r"""Plot lines of constant mixing ratio.
#
#        Adds lines of constant mixing ratio (isohumes) to the
#        plot. The default style of these lines is dashed green lines with an
#        alpha value of 0.8. These can be overridden using keyword arguments.
#
#        Parameters
#        ----------
#        w : array_like, optional
#            Unitless mixing ratio values to plot. If none are given, default
#            values are used.
#        p : array_like, optional
#            Pressure values to be included in the isohumes. If not
#            specified, they will be linearly distributed across the current
#            plotted pressure range up to 600 mb.
#        kwargs
#            Other keyword arguments to pass to :class:`matplotlib.collections.LineCollection`
#
#        Returns
#        -------
#        matplotlib.collections.LineCollection
#            instance created
#
#        See Also
#        --------
#        :class:`matplotlib.collections.LineCollection`
#
#        """
#        # Default mixing level values if necessary
#        if w is None:
#            w = np.array([0.0004, 0.001, 0.002, 0.004, 0.007, 0.01,
#                          0.016, 0.024, 0.032]).reshape(-1, 1)
#
#        # Set pressure range if necessary
#        if p is None:
#            p = np.linspace(600, max(self.ax.get_ylim())) * units.mbar
#
#        # dewpoint lines, drawn on Celcius axis unless the xaxis has other units specified
#        # or if the xlimits suggest a kelvin range (heuristically)
#        td = dewpoint(vapor_pressure(p, w))
#        if self.ax.xaxis.units is not None:
#            td=td.to(self.ax.xaxis.units)
#        elif self.ax.get_xlim()[0] > 150:
#            td=td.to(units.degK)
#        
#        # Assemble data for plotting
#        linedata = [np.vstack((t, p)).T for t in td]
#
#        # Add to plot
#        kwargs.setdefault('colors', 'g')
#        kwargs.setdefault('linestyles', 'dashed')
#        kwargs.setdefault('alpha', 0.8)
#        return self.ax.add_collection(LineCollection(linedata, **kwargs))