# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 11:47:15 2019
    Examine AWS dataset(s)
@author: jgreensl
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas
import iris
from datetime import datetime, timedelta

from utilities import fio, plotting, utils

_sn_='AWS'

# Wagerup AWS entries to plot
# there may be a better way to set all these plotting parameters...
__AWS_PLOTS__ = {'Wind direction':{'dfnames':['Dta10','Dta30'],
                                   'colors':['k','darkgrey'],
                                   'linestyles':['None','None'],
                                   'markers':['.','^']},
                 'Wind speed':{'dfnames':['S10','S30'],
                               'colors':['k','darkgrey'],
                               'linestyles':['None','None'],
                               'markers':['.','^']},
                 #'Gusts':['SX10','SX30'],
                 'Pressure':{'dfnames':['QFE'],
                             'colors':['k',],
                             'linestyles':['None',],
                             'markers':['.']},
                 'Relative humidity':{'dfnames':['RH'],
                                      'colors':['k'],
                                      'linestyles':['None'],
                                      'markers':['.']},
                 'Solar radiation':{'dfnames':['SR'],
                                    'colors':['k'],
                                    'linestyles':['None'],
                                    'markers':['.']},
                 'Temperature':{'dfnames':['T','T10','T30'],
                                'colors':['k','darkgrey','grey'],
                                'linestyles':['None','None','None'],
                                'markers':['v','.','^']},
                }

def df_time_series(df, subplots=None, units=None, legend=True):
    """
    plot time series from df input
    subplots = dict: {'title1':['columna','columnb',...],...}
    if subplots is defined, subplot i will show df entries listed by subplots[i]
    """
    if subplots is None:
        subplots = { df.keys()[i]:df.keys()[i] for i in range(len(df.keys())) }
    
    f, axes = plt.subplots(nrows=len(subplots), figsize=(12, 24))
    
    for i, title in enumerate(subplots.keys()):
        ax=axes[i]
        plt.sca(ax)
        #print("DEBUG:", title, subplots[title])
        dfd = subplots[title]
        
        for ii, name in enumerate(dfd['dfnames']):
            #print("DEBUG:", name, ii, dfd['markers'], dfd['linestyles'], dfd['colors'])
            df[name].plot(alpha=0.5, ax=ax, 
                          marker=dfd['markers'][ii],
                          linestyle=dfd['linestyles'][ii],
                          color=dfd['colors'][ii])
        
        if legend:
            plt.legend()
        
        #df[dfd['dfnames']].plot(marker=dfd['markers'], alpha=0.5, 
        #                        linestyles=dfd['linestyles'], ax=ax,
        #                        colors=dfd['colors'])
        plt.grid(which='major',axis='x', alpha=0.6)
        if ax != axes[-1]:
            plt.xlabel('')
            ax.tick_params(axis='x',labelbottom='off')
        else:
            #plt.xlabel('time (AWST)')
            plt.xlabel(df.index.name)
        plt.title(title,fontsize=15)
        if units is not None:
            plt.ylabel(units[title],fontsize=12)
    return f, axes



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
    units = {key:dfa[subplots[key]['dfnames'][0]]['units'] for key in subplots.keys()}

    df_time_series(dfsub,subplots=subplots,units=units)
    plt.suptitle('Wagerup AWS '+substr,fontsize=19)
    fio.save_fig('measurements',_sn_,'summary_wagerup',plt)

def compare_site_to_model(AWS='wagerup',
                          model_run='waroona_run1',
                          dtoffset=timedelta(hours=8), # AWST offset
                          heights=[1,10,30], # what heights to look at in model output
                          umhours=None, # load model data for these hours
                          showrange=None, # subset the time series to this range
                          ):
    """
    AWS='wagerup', model_run='waroona_run1',dtoffset=timedelta(hours=8)):
    
    compare site to model run
    """
    
    lat,lon = plotting._latlons_['wagerup']
    extent = [lon-.02, lon+.02, lat-.02, lat+.02] # just grab real close to latlon

    ## Read AWS: 
    if AWS == 'wagerup':
        df_aws, aws_attrs = fio.read_AWS_wagerup() 

    # Default model hours is ALL
    if umhours is None:
        umhours = fio.model_outputs[model_run]['filedates']

    ## Read Model output:
    data_model0 = fio.read_model_run(model_run, fdtime=umhours, extent=extent, add_topog=True, add_winds=True, add_RH=True, add_z=True)
    
    wanted_cube_names=['wind_direction','s','air_pressure','relative_humidity','air_temperature', 'z_th', 'surface_altitude']
    wanted_cube_units=[None,            'm s-1', 'hPa',        '%',             'deg_c' ,       None,     None]
    #Just want time and level at our location
    data_model = data_model0.extract(wanted_cube_names)
    for i in range(len(data_model)):
        data_model[i] = data_model[i].interpolate([('longitude',lon), ('latitude',lat)],
                                                  iris.analysis.Linear())
        if wanted_cube_units[i] is not None:
            data_model[i].convert_units(wanted_cube_units[i])
    
    # which model levels represent 10 and 30m?
    z, topog = data_model.extract(['z_th','surface_altitude'])
    # Just keep time series in our cubelist
    data_model.remove(z)
    data_model.remove(topog)
    
    # height above ground level
    agl = z.data - topog.data
    agl_coord = iris.coords.AuxCoord(agl, var_name='height_agl', units='m')
    agl_coord.guess_bounds()
    
    # interpolate onto 1m, 10m, and 30m levels
    dict_model = {}
    for i in range(len(data_model)):
        data_model[i].add_aux_coord(agl_coord,1)
        # interpolate to 1, 10, 30
        model_cube = data_model[i].interpolate([('height_agl',heights)], iris.analysis.Linear())
        
        # make time series for each vertical level in a dictionary
        for zi in range(len(heights)):
            cubedata = model_cube[:,zi].data
            if isinstance(cubedata, np.ma.MaskedArray):
                cubedata = cubedata.data
            dict_model[model_cube.name()+str(heights[zi])]=cubedata
        
    
    dtstr = [(dt+dtoffset).strftime("%Y-%m-%d %H:%M:%S") for dt in utils.dates_from_iris(data_model[0])]
    dict_model['WAST']=dtstr
        
    df_model = pandas.DataFrame(dict_model)
    df_model['WAST'] = pandas.to_datetime(df_model['WAST'])
    # set the datetime column to be the index
    df_model = df_model.set_index('WAST')
    
    # Merge the cubes
    # outer join gets the union, inner join gets intersection
    df_both = pandas.merge(df_aws,df_model, how='outer', left_index=True, right_index=True)
    
    df_both_subset = df_both
    if showrange is not None:
        rangestr = [showrange[0].strftime("%Y-%m-%d %H"),
                     showrange[1].strftime("%Y-%m-%d %H")]
        df_both_subset = df_both.loc[rangestr[0]:rangestr[1]]
    
    ## Update subplots dictionnary with model column names
    subplots = __AWS_PLOTS__
    for title, mname in zip(['Wind direction','Wind speed','Pressure','Relative humidity','Temperature'],
                            ['wind_direction','s','air_pressure','relative_humidity','air_temperature']):
        # add model name to list of dataframe columns to be plotted within corresponding subplot
        subplots[title]['dfnames'].extend([mname+'%d'%d for d in [1,10,30]])
        # add marker,color,linestyle too
        subplots[title]['colors'].extend(['red','magenta','darkpink'])
        subplots[title]['linestyles'].extend(['None']*3)
        subplots[title]['markers'].extend(['v','.','^'])
    
    ## now call method which plots the dataframe directed by my dictionary of subplot names
    df_time_series(df_both_subset, subplots=subplots)
    fio.save_fig(model_run=model_run, script_name=_sn_, 
                 pname='%s_vs_%s'%(AWS,model_run), plt=plt)

if __name__=='__main__':
    summary_wagerup()
    compare_site_to_model()