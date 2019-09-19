# -*- coding: utf-8 -*-

def ax_skewt(tlims=[240,330],plims=[1050,100], th_bl=None, q_bl=None):
    '''
    Using Jeff's saturationpoint code, create base for log p skew t plot
    EG: 
        ax=ax_skewt(tlims=[250,375],plims=[1050,50])
        ax.semilogy(T[0,:,50,50].data, p[0,:,50,50].data/100., 'k')
        plt.show()
    Inputs: 
        temperature limits, pressure limits for plotting
    Returns:
        axis with stuff already drawn
    Additional: 
        boundary layer theta and mixing ratio can be added to draw SP curve 
        
    
    '''
    
    #import matplotlib as mpl
    #import matplotlib.pyplot as plt
    #from matplotlib.ticker import ScalarFormatter#, MultipleLocator
    #from matplotlib.projections import register_projection
    #import numpy as np
    
    #import skewt
    #import thermo
    
    # Dry adiabats rise from right to left, representing constant potential temperature
    def plot_dryadiabat(ax,theta,p,p0):
        # p should be an array of pressures at which the adiabat is plotted
        T = theta*(p/p0)**thermo.kappa
        ax.semilogy(T, p, 'g-',alpha=0.5)
        
    # mixing ratios rise left to right
    def plot_mixrat(ax,r,p):
        # p should be an array of pressures (hPa) at which the adiabat is plotted
        e = p * r / (thermo.eps + r)
        loge = np.log(e)
        Td = (243.5*loge - 440.8)/(19.48 - loge) + 273.15
        ax.semilogy(Td,p,'g-',alpha=0.5)
    
    def plot_moist_adiabats(ax, t0=None, p=None, **kwargs):
        r"""Plot moist adiabats.
        
        Temperatures in K, pressures in hPa.

        Adds saturated pseudo-adiabats (lines of constant equivalent potential
        temperature) to the plot. The default style of these lines is dashed
        blue lines with an alpha value of 0.5. These can be overridden using
        keyword arguments.

        Parameters
        ----------
        t0 : array_like, optional
            Starting temperature values in Kelvin. If none are given, they will be
            generated using the current temperature range at the bottom of
            the plot.
        p : array_like, optional
            Pressure values to be included in the moist adiabats. If not
            specified, they will be linearly distributed across the current
            plotted pressure range.
        kwargs
            Other keyword arguments to pass to :class:`matplotlib.collections.LineCollection`

        Returns
        -------
        matplotlib.collections.LineCollection
            instance created

        See Also
        --------
        :func:`~metpy.calc.thermo.moist_lapse`
        :meth:`plot_dry_adiabats`
        :class:`matplotlib.collections.LineCollection`

        """
        # Determine set of starting temps if necessary
        if t0 is None:
            xmin, xmax = ax.get_xlim()
            t0 = np.arange(240, xmax, 10)
        print("DEBUG: t0",t0.shape, t0)
        # Get pressure levels based on ylims if necessary
        if p is None:
            p = np.linspace(1050,100, 96)
        print("DEBUG: p",p.shape, p)
        # Assemble into data for plotting
        t = thermo.moist_lapse(p*100, t0[:, np.newaxis], 100000) # in Kelvin
        linedata = [np.vstack((ti, p)).T for ti in t]
        print("DEBUG: t",t.shape, t)
        # Add to plot
        kwargs.setdefault('colors', 'b')
        kwargs.setdefault('linestyles', 'dashed')
        kwargs.setdefault('alpha', 0.5)
        #return ax.add_collection(LineCollection(linedata, **kwargs))
        ax.add_collection(LineCollection(linedata, **kwargs))
    
    #    def plot_moistadiabat(ax, T, p, p0=1e5):
    #        # moist adiabat, 
    #        #at the saturated mixing ratios (kg/kg) over pressure (Pa)
    #        
    #        #for Tk in np.array(T):
    #        r_sat = thermo.r_sat(T,p0)# kg/kg = f(K,Pa)
    #        th_e = thermo.theta_e(T, r_sat, p0) # K/Pa = f(K, kg/kg, Pa)
    #        print("DEBUG: th_e",th_e.shape, th_e)
    #        print("DEBUG: p",p.shape, p)
    #        print("DEBUG: r_sat",r_sat.shape,r_sat)
    #        e = p * r_sat / (thermo.eps + r_sat) 
    #        loge = np.log( 0.01*np.maximum(e, 1e-100) )
    #        T_lcl = 2840.0 / (3.5*np.log(T[0]) - loge - 4.805) + 55.0
    #        #print("DEBUG:",th_e.shape, p.shape, r_sat.shape, T_lcl.shape)
    #        Te = ( th_e * (100000/p0) ** (.079912*r_sat - thermo.kappa) 
    #              * np.exp((2.54-3376/T_lcl)*r_sat*(1+.81*r_sat)))
    #        print("DEBUG: Te",Te.shape, Te)
    #        ax.semilogy(Te,p/100.,'g:',linewidth=2)
    #        #for Tk in np.array(T):
    #        #    r_sat = thermo.r_sat(Tk,p) # kg/kg = f(K,Pa)
    #        #    
    #        #    th_e = thermo.theta_e(Tk, r_sat, p) # K/Pa = f(K, kg/kg, Pa)
    #        #    e = p * r_sat / (thermo.eps + r_sat) 
    #        #    loge = np.log( 0.01*np.maximum(e, 1e-100) )
    #        #    T_lcl = 2840.0 / (3.5*np.log(Tk) - loge - 4.805) + 55.0
    #        #    #print("DEBUG:",th_e.shape, p.shape, r_sat.shape, T_lcl.shape)
    #        #    Te = ( th_e * (100000/p) ** (.079912*r_sat - thermo.kappa) 
    #        #          * np.exp((2.54-3376/T_lcl)*r_sat*(1+.81*r_sat)))
    #        #    ax.semilogy(Te,p,'g:')
    #        #print("DEBUG: Te",Te)
    #        #print("DEBUG: p",p)
    #        #print("DEBUG: r_sat",r_sat)
    #        #print("DEBUG: T_lcl",T_lcl)
        
    matplotlib.rcParams['font.size'] = 18.0
    matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
    mtfs = 22.0  # fontsize for math text
    
    register_projection(skewt.SkewXAxes)
    
    p0 = 1e5
    
    fig = plt.figure(1)#figsize=(6.5875, 6.2125))
    ax = fig.add_subplot(111, projection='skewx')
    
    yticks = np.ravel(np.outer(10**np.arange(0,3),np.linspace(1,10,10)))
    
    plt.grid(True)
    for thx in range(300,2000,100):
        plot_dryadiabat(ax,thx,yticks,1e3)
    for rx in [1e-4,2e-4,4e-4,1e-3,2e-3,4e-3,1e-2,2e-2,4e-2,0.1,0.2,0.4,1.0]:
        plot_mixrat(ax,rx,yticks)
    plot_moist_adiabats(ax)
    #for Tx in range(270,325,10):
    #    Tarr = np.zeros([20])+Tx
    #    parr = np.logspace(5,4,20)
    #    plot_moistadiabat(ax, Tarr,parr)
    
    
    
    if th_bl is not None and q_bl is not None:
        #th_bl = 303
        #q_bl = 0.005
        
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
    
        # Plot the data using normal plotting functions, in this case using
        # log scaling in Y, as dicatated by the typical meteorological plot
        ax.semilogy(T_sp, p_sp*1e-2, 'r')
        #ax.semilogy([T_fire,T_sp[-1],Td_fire],[p0*1e-2,p_sp[-1]*1e-2,p0*1e-2],'b-')
        ax.semilogy([T_bl,T_sp[0],Td_bl],[p0*1e-2,p_sp[0]*1e-2,p0*1e-2],'b-')
        
        ax.text(0.8,0.9,  r'$\theta_{bl}='  +'{:.0f}'.format(th_bl)   +'\mathrm{K}$',transform=ax.transAxes,fontsize=mtfs)
        ax.text(0.8,0.825,r'$q_{bl}='       +'{:.1f}'.format(q_bl*1e3)+'\mathrm{g/kg}$',transform=ax.transAxes,fontsize=mtfs)
        ax.text(0.8,0.75, r'$\gamma='       +'{:.1f}$'.format(gamma),transform=ax.transAxes,fontsize=mtfs)
        ax.text(0.8,0.675,r'$\delta='       +'{:.1f}$'.format(delta*1e-3),transform=ax.transAxes,fontsize=mtfs)
        ax.text(0.8,0.6,  r'$\theta_{fire}='+'{:.0f}'.format(th_fire)+'\mathrm{K}$',transform=ax.transAxes,fontsize=mtfs)
        ax.text(0.8,0.525,r'$q_{fire}='     +'{:.1f}'.format(q_fire*1e3)+'\mathrm{g/kg}$',transform=ax.transAxes,fontsize=mtfs)
    
    # Disables the log-formatting that comes with semilogy
    ax.yaxis.set_major_formatter(tick.ScalarFormatter())
    ax.set_yticks(yticks)
    ax.set_ylim(plims[0], plims[1])
    
    ax.xaxis.set_major_locator(tick.MultipleLocator(10))
    #ax.set_xlim(Td_bl,400)
    ax.set_xlim(tlims[0],tlims[1])
    
    ax.set_xlabel('$T \, (\mathrm{K})$')
    ax.set_ylabel('$p \, (\mathrm{hPa})$')
    
    fig.set_size_inches(20,10)
    #plt.show()
    #plt.savefig('skewT_th{:.0f}_q{:.0f}_gam{:.0f}_del{:.0f}.png'.format(th_bl,q_bl*1e3,gamma,delta*1e-3),dpi=200)
    return ax 