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
from copy import deepcopy
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
def summary_wagerup(d0=None,dN=None, UTC=True):
    """
    show wagerup AWS time series
    optionally subset time to day: d0, or from d0 to dN inclusive
    """
    df, dfa = fio.read_AWS_wagerup(UTC)
    
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
    subplots = deepcopy(__AWS_PLOTS__)
    units = {key:dfa[subplots[key]['dfnames'][0]]['units'] for key in subplots.keys()}

    df_time_series(dfsub,subplots=subplots,units=units)
    plt.suptitle('Wagerup AWS '+substr,fontsize=19)
    fio.save_fig('measurements',_sn_,'summary_wagerup'+['','_UTC'][UTC],plt)

def combine_site_and_model(AWS='wagerup', model_run='waroona_run1',
                           heights=[1,10,30],
                           groupbystr=None,
                           UTC=True):
    """
    read AWS dataframe, add model columns to dataframe
        optionally resample intervals (mean) using groupbystr
    
    arguments
    ---------
    heights: list, interpolate model parameters at site to these heights (m)
    groupbystr: string, make all columns binned (by '30Min' for example)
    """
    lat,lon = plotting._latlons_[AWS]
    extent = [lon-.02, lon+.02, lat-.02, lat+.02] # just grab real close to latlon

    ## Read AWS: eventually handle multiple weather stations
    #if AWS == 'wagerup':
    subplots = deepcopy(__AWS_PLOTS__)
    df_aws, aws_attrs = fio.read_AWS_wagerup(UTC=UTC)

    umhours = fio.model_outputs[model_run]['filedates']
    
    ## Read Model output:
    data_model0 = fio.read_model_run(model_run, fdtime=umhours, extent=extent, 
                                     add_topog=True, add_winds=True, 
                                     add_RH=True, add_z=True)
    #print("DEBUG:",data_model0)
    wanted_cube_names=['wind_direction', 's', 'air_pressure', 
                       'relative_humidity', 'air_temperature', 'z_th', 
                       'surface_altitude']
    wanted_cube_units={'s':'m s-1',
                       'air_pressure':'hPa',
                       'relative_humidity':'%',
                       'air_temperature':'deg_c'
                       }
    #Just want time and level at our location
    data_model = data_model0.extract(wanted_cube_names,strict=True) # strict means one cube per constraint
    for i in range(len(data_model)):
        data_model[i] = data_model[i].interpolate([('longitude',lon), ('latitude',lat)],
                                                  iris.analysis.Linear())
        dname = data_model[i].name()
        if dname in wanted_cube_units.keys():
            #print("DEBUG: converting",dname, "from ",data_model[i].units," to ",wanted_cube_units[dname])
            data_model[i].convert_units(wanted_cube_units[dname])
    
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
        
        # put time series for each vertical level in a dictionary
        for zi in range(len(heights)):
            cubedata = model_cube[:,zi].data
            if isinstance(cubedata, np.ma.MaskedArray):
                cubedata = cubedata.data
            dict_model[model_cube.name()+str(heights[zi])]=cubedata
        
    dt_model0 = utils.dates_from_iris(data_model[0]) # UTC
    dt_model=dt_model0
    if not UTC:
        dt_model = [ dt + timedelta(hours=8) for dt in dt_model0]
        
    dtstr = [(dt).strftime("%Y-%m-%d %H:%M:%S") for dt in dt_model]
    dtname = ['WAST','UTC'][UTC]
    dict_model[dtname]=dtstr
        
    df_model = pandas.DataFrame(dict_model)
    df_model[dtname] = pandas.to_datetime(df_model[dtname])
    # set the datetime column to be the index
    df_model = df_model.set_index(dtname)
    
    # Merge the cubes
    # outer join gets the union, inner join gets intersection
    df_both = pandas.merge(df_aws,df_model, how='outer', left_index=True, right_index=True)
    
    ## Update subplots dictionnary with model column names
    
    for title, mname in zip(['Wind direction','Wind speed','Pressure','Relative humidity','Temperature'],
                            ['wind_direction','s','air_pressure','relative_humidity','air_temperature']):
        # add model name to list of dataframe columns to be plotted within corresponding subplot
        subplots[title]['dfnames'].extend([mname+'%d'%d for d in [1,10,30]])
        # add marker,color,linestyle too
        subplots[title]['colors'].extend(['crimson','fuchsia','blueviolet'])
        subplots[title]['linestyles'].extend(['None']*3)
        subplots[title]['markers'].extend(['v','.','^'])
    
    df = df_both
    if groupbystr is not None:
        df = df_both.resample('30Min')
    
    return df, subplots
    
def compare_site_to_model(AWS='wagerup',
                          model_run='waroona_run1',
                          heights=[1,10,30], # what heights to look at in model output
                          showrange=None, # subset the time series to this range
                          UTC=True,
                          ):
    """
    AWS='wagerup', model_run='waroona_run1',dtoffset=timedelta(hours=8)):
    
    compare site to model run
    """
    df, subplots = combine_site_and_model(AWS=AWS, model_run=model_run,
                                          heights=heights, UTC=UTC)
    
    df_subset = df
    if showrange is not None:
        rangestr = [showrange[0].strftime("%Y-%m-%d %H"),
                     showrange[1].strftime("%Y-%m-%d %H")]
        df_subset = df.loc[rangestr[0]:rangestr[1]]
    
    print("DEBUG:",df_subset)
    ## now call method which plots the dataframe directed by my dictionary of subplot names
    df_time_series(df_subset, subplots=subplots)
    fio.save_fig(model_run=model_run, script_name=_sn_, 
                 pname='%s_%s_all%s'%(AWS,model_run,['','_UTC'][UTC]), plt=plt)
    
    ## aggregate and compare model to aws
    df_agg = df.resample('30Min').mean().dropna()
    print("DEBUG:",df_agg)
    
    df_time_series(df_agg, subplots=subplots)
    fio.save_fig(model_run=model_run, script_name=_sn_, 
                 pname='%s_%s_aggregate%s'%(AWS,model_run,['','_UTC'][UTC]), plt=plt)
    
    
    
    return df_agg, subplots
    
if __name__=='__main__':
    d0,dN = datetime(2016,1,4,10), datetime(2016,1,7,10)
    summary_wagerup(d0,dN,UTC=True)
    summary_wagerup(d0,dN,UTC=False)
    
    for mr in ['waroona_run1','waroona_old']:
        for UTC in [True, False]:
            compare_site_to_model(AWS='wagerup',
                                  model_run=mr,
                                  showrange=[datetime(2016,1,5,6),datetime(2016,1,7,1)],
                                  UTC=UTC)

#(df1,sp1), (dfold,spold) = df

