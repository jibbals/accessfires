#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker, colors, patches
import numpy as np
from netCDF4 import Dataset
import pandas

from datetime import datetime, timedelta
import timeit
import warnings

import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.io.img_tiles as cimgt

from utilities import utils,fio,plotting


# Try reading file using iris
import iris
import iris.quickplot as qplt
from iris.experimental.equalise_cubes import equalise_attributes


# Wagerup AWS entries to plot
__AWS_PLOTS__ = {'Wind direction':['Dta10','Dta30'],
                 'Wind speed':['S10','S30'],
                 #'Gusts':['SX10','SX30'],
                 'Pressure':['QFE'],
                 'Relative humidity':['RH'],
                 'Solar radiation':['SR'],
                 'Temperature':['T','T10','T30'],
                }

### CODE


def df_time_series(df, subplots=None, marker='.', linestyle='None', units=None):
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
        print("DEBUG:",title, subplots[title])
        df[subplots[title]].plot(marker=marker, alpha=0.5, linestyle=linestyle, ax=ax)
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

"""
AWS='wagerup', model_run='waroona_run1',dtoffset=timedelta(hours=8)):

compare site to model run
"""
AWS='wagerup'
model_run='waroona_run1'
dtoffset=timedelta(hours=8) # AWST offset
lat,lon = plotting._latlons_['wagerup']
extent = [lon-.02, lon+.02, lat-.02, lat+.02] # just grab real close to latlon
compare_list = ['Wind direction','Wind speed', 'Pressure', 'Relative humidity','Temperature']

## Read AWS: 
if AWS == 'wagerup':
    df_aws, aws_attrs = fio.read_AWS_wagerup() 

## Read Model output: need wind direction, relative humidity to be added
data_model0 = fio.read_model_run(model_run, fdtime=datetime(2016,1,5,15), extent=extent, add_topog=True, add_winds=True, add_RH=True, add_z=True)
print(data_model0)

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
data_model.remove(z)
data_model.remove(topog)
print("DEBUG:", data_model)

# height above ground level
agl = z.data - topog.data
agl_coord = iris.coords.AuxCoord(agl, var_name='height_agl', units='m')
agl_coord.guess_bounds()

# interpolate onto 1m, 10m, and 30m levels
dict_model = {}
heights = [1,10,30]
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

print(dict_model)

df_model = pandas.DataFrame(dict_model)
df_model['WAST'] = pandas.to_datetime(df_model['WAST'])
# set the datetime column to be the index
df_model = df_model.set_index('WAST')

# Merge the cubes?
# outer join gets the union, inner join gets intersection
df_both = pandas.merge(df_aws,df_model, how='outer', left_index=True, right_index=True)
print(df_both.loc['2016-01-05 23'])
df_both_subset = df_both.loc['2016-01-05 18':'2016-01-06 05']

#for cube,cubename in zip([wdc,wsc,pc,rhc,tc],['])
#wd = wdc.

#'Gusts':['SX10','SX30'],
#'Pressure':['QFE'],
#'Relative humidity':['RH'],
#'Solar radiation':['SR'],
#'Temperature':['T','T10','T30'],

__AWS_PLOTS__['Wind direction'].extend(['wind_direction%d'%d for d in [1,10,30]])
__AWS_PLOTS__['Wind speed'].extend(['s%d'%d for d in [1,10,30]])
__AWS_PLOTS__['Pressure'].extend(['air_pressure%d'%d for d in [1,10,30]])
__AWS_PLOTS__['Relative humidity'].extend(['relative_humidity%d'%d for d in [1,10,30]])
__AWS_PLOTS__['Temperature'].extend(['air_temperature%d'%d for d in [1,10,30]])

df_time_series(df_both_subset, subplots=__AWS_PLOTS__)
# Create dataframe with date indices

# merge with AWS dataframe
# Send to the plotting hole