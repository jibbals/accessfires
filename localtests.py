#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker, colors, patches, dates
import numpy as np


from datetime import datetime, timedelta

from utilities import fio, plotting, utils

_sn_ = 'localtests'



import pandas
import AWS_series
"""
AWS='wagerup',
                      model_run='waroona_run1',
                      heights=[1,10,30], # what heights to look at in model output
                      showrange=None, # subset the time series to this range
                      UTC=True,
                      ):
AWS='wagerup', model_run='waroona_run1',dtoffset=timedelta(hours=8)):

compare site to model run
"""
AWS='wagerup'
model_run='waroona_run1'
heights=[1,10,30] # what heights to look at in model output
showrange=None # subset the time series to this range
UTC=True


df, subplots = AWS_series.combine_site_and_model(AWS=AWS, model_run=model_run,
                                      heights=heights, UTC=UTC)

df_subset = df
if showrange is not None:
    rangestr = [showrange[0].strftime("%Y-%m-%d %H"),
                 showrange[1].strftime("%Y-%m-%d %H")]
    df_subset = df.loc[rangestr[0]:rangestr[1]]

## now call method which plots the dataframe directed by my dictionary of subplot names
AWS_series.df_time_series(df_subset, subplots=subplots)
fio.save_fig(model_run=model_run, script_name=_sn_, 
             pname='%s_%s_all%s'%(AWS,model_run,['','_UTC'][UTC]), plt=plt)

## aggregate and compare model to aws
df_agg = df.resample('30Min').mean()
AWS_series.df_time_series(df_agg, subplots=subplots)
fio.save_fig(model_run=model_run, script_name=_sn_, 
             pname='%s_%s_aggregate%s'%(AWS,model_run,['','_UTC'][UTC]), plt=plt) 

## Now we look at statistics
# if any columns are NAN, drop whole row
df_agg = df_agg.dropna(how='any') 

# diurnal half hourly cycle
diurnal = df.resample('30Min').mean() # 30 min averages
diurnal = diurnal.groupby([diurnal.index.hour, diurnal.index.minute]).mean() # combine multiple days

# now get correlation,bias between certaint column pairs:
titles = ['Wind dir', 'RH','T','Wind speed','Pressure']
pairs = [['Dta10','wind_direction10'],['RH','relative_humidity10'],['T','air_temperature1'],['S10','s10'],['QFE','air_pressure10']]

plt.close()
fig, axes = plt.subplots(nrows=5, sharex=True, figsize=[7,12])

# for each pair, show them detrended and the correlation before/after detrending
for j, (pair, title) in enumerate(zip(pairs,titles)):
    plt.sca(axes[j])
    
    corr = df_agg[pair].corr()
    
    # detrend the pair:
    # model data being detrended with measurements (longer time series)
    dpair = [pair[0],pair[0]] # just measured column twice
    hours,minutes = df_agg.index.hour, df_agg.index.minute
    diff = []
    for jj in range(len(df_agg)):
        H,M = hours[jj], minutes[jj]
        diff.append( df_agg[pair].values[jj] - diurnal[dpair].loc[(H,M)].values)
    # make dataframe using differences, column names, and set index to match df_agg
    detrended = pandas.DataFrame(diff, columns=pair)
    detrended.index = df_agg.index
    # correlation
    detrended_corr = detrended.corr()
    
    #detrended.plot(colors=['k','m'])
    line1 = plt.plot_date(detrended.index, detrended[pair[0]].values, 
                          linestyle='-', color='k', marker='None')
    line2 = plt.plot_date(detrended.index, detrended[pair[1]].values, 
                          linestyle='-', color='m', marker='None')
    
    # add legend with two extra bits of info in it...
    handles = [line1[0],line2[0]]
    labels = [pair[0], pair[1]]
    # add handle,label with correlations in it
    handles.append(patches.Patch(color='white'))
    labels.append('r = %.3f'%corr.values[0,1])
    handles.append(patches.Patch(color='white'))
    labels.append('dt-r = %.3f'%detrended_corr.values[0,1])
    plt.legend(handles,labels,loc='best')
    
    # correlations before and after detrending
    #plt.annotate("trendy    r = %.3f"%corr.values[0,1], 
    #             xy=(.085,.75), xycoords='axes fraction')
    #plt.annotate("detrended r = %.3f"%detrended_corr.values[0,1], 
    #             xy=(.085,.7), xycoords='axes fraction')
    
    if j < 4:
        # remove xticks etc
        plt.xticks([],[])
        plt.xlabel('')
        plt.gca().axes.xaxis.set_visible(False)
plt.suptitle('Detrended comparable metrics',fontsize=20)
plt.subplots_adjust(top=.95,hspace=0.02)
# rotate and align the tick labels so they look better
fig.autofmt_xdate()
plt.show()