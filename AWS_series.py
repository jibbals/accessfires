# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 11:47:15 2019
    Examine AWS dataset(s)
@author: jgreensl
"""

import numpy as np
import matplotlib.pyplot as plt

from utilities import fio

_sn_='AWS'




# read wagerup
def summary_wagerup(d0=None,dN=None):
    """
    show wagerup AWS time series
    optionally subset time to day: d0, or from d0 to dN inclusive
    """
    df, dfa = fio.read_AWS_wagerup()
    
    dfsub=df
    substr=''
    if d0 is not None and dN is not None:
        substr=d0.strftime("%Y-%m-%d %H:%M")+' - '+dN.strftime("%Y-%m-%d %H:%M")
        dfsub=df.loc[d0.strftime("%Y-%m-%d %H:%M"):dN.strftime("%Y-%m-%d %H:%M")]
    elif d0 is not None:
        # Just pull out one day
        substr = d0.strftime("%Y-%m-%d %H:%M")
        dfsub = df.loc[d0.strftime("%Y-%m-%d %H:%M")]
    
    # plot some stuff
    subplots = __AWS_PLOTS__
    plt.close()
    f, axes = plt.subplots(nrows=len(subplots), figsize=(12, 24))
    
    for i, title in enumerate(subplots.keys()):
        ax=axes[i]
        plt.sca(ax)
        dfsub[subplots[title]].plot(marker='.', alpha=0.5, linestyle='None', ax=ax)
        plt.grid(which='major',axis='x', alpha=0.6)
        if ax != axes[-1]:
            plt.xlabel('')
            ax.tick_params(axis='x',labelbottom='off')
        else:
            plt.xlabel('time (AWST)')
        plt.title(title,fontsize=15)
        plt.ylabel(dfa[subplots[title][0]]['units'],fontsize=12)
    plt.suptitle('Wagerup AWS '+substr,fontsize=19)
    fio.save_fig('measurements',_sn_,'summary_wagerup',plt)

def compare_site(AWS='wagerup', model_run='waroona_run1'):
    """
    compare site to model run
    """
    lat,lon = plotting._latlons_['wagerup']
    extent = [lon-.02, lon+.02, lat-.02, lat+.02] # just grab real close to latlon
    compare_list = ['Wind direction','Wind speed', 'Pressure', 'Relative humidity','Temperature']
    
    ## Read AWS: 
    if AWS == 'wagerup':
        data_aws, aws_attrs = fio.read_AWS_wagerup() 
    
    ## Read Model output: need wind direction, relative humidity to be added
    data_model = fio.read_model_run(model_run, extent=extent, add_winds=True, add_RH=True)
    # Pull out time series for Wind dir, Wind speed, Pressure, RH, Temp
    wd,ws,p,rh,t = data_model.extract(['wind_direction','s','air_pressure','relative_humidity','air_temperature'])
    # Create dataframe with date indices
    # merge with AWS dataframe
    # Send to the plotting hole
    

if __name__=='__main__':
    summary_wagerup()