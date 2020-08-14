#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:55:32 2019

  Script to run tests etc before locating them properly
  
@author: jesse
"""


# IMPORTS

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches, colors, cm
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection
from datetime import datetime, timedelta
import iris # so we can add level number constraint
import iris.analysis
import iris.quickplot as qplt

# read tiff file
from osgeo import gdal, osr
import cartopy.crs as ccrs

from utilities import utils, plotting, fio, constants

from scipy import interpolate

import pandas

# Dict to hold all the AWS metadata
_AWS_={
    "moree_airport":{
        "path_AIFS":"data/AWS/MoreeAirport.csv",
        "latlon":[-29.4946, 149.8505], # degrees
        "altitude":214, # metres ASL
        },
    "coonabarabran_airport":{
        "path_AIFS":"data/AWS/CoonabarabranAirport.csv",
        "latlon":[-31.3330,149.2699],
        "altitude":645,
        },
    }
## Read AIFS based csv files
def AIFS_read_path(path):
    """
    Return path into Pandas table:
        Station	DF	Cur %	Grass t/ha	EST HH:MM dd/mm/yy	T	Td	RH	Wd	Ws km/h	Wg km/h	FFDR/FFDI	GFDR/GFDI

    """
    
    ## Read data
    print("INFO: reading ",path)
    data = pandas.read_csv(path) # skiprows=1 removes titles

    # read datetimes
    dind_str='EST HH:MM dd/mm/yy' # local time column
    data[dind_str] = pandas.to_datetime(data[dind_str], dayfirst=True)
    
    # Wind direction from string to degrees clockwise from due North
    wind_dir_map = {"N":0,"NNE":22.5,"NE":45,"ENE":67.5,"E":90,
                    "ESE":112.5,"SE":135,"SSE":157.5,"S":180,
                    "SSW":202.5,"SW":225,"WSW":247.5,"W":270,
                    "WNW":292.5,"NW":315,"NNW":337.5,}
    data['Wd'] = data['Wd'].map(wind_dir_map)
    
    def number_from_fdi(fdistr): 
        return int(fdistr.split(' ')[1])
    for key in ['FFDR/FFDI','GFDR/GFDI']:
        data[key] = data[key].apply(number_from_fdi,convert_dtype=True)
    
    
    # make sure data columns are numeric (not string)
    for key in data.keys():
        if key not in ['Station',dind_str]:
            data[key] = pandas.to_numeric(data[key], errors='coerce')
    
    # set local time to be the index, and sort by date
    data = data.set_index(dind_str)
    data = data.sort_values(by=dind_str)
    return data

def AWS_plot_timeseries(df,key,d0=None,dN=None,**plotargs):
    """
    ARGS:
        df: dataframe with datetime index
        key: string name of item from df to be plotting
    """
    
    subdf=df.copy()
    if d0 is not None:
        subdf = subdf.loc[d0:dN]
    
    subdf[key].plot(**plotargs)

## READ MOREE
sname='moree_airport'
path = _AWS_[sname]['path_AIFS']
latlon = _AWS_[sname]['latlon']

aws=AIFS_read_path('data/AWS/MoreeAirport.csv')
print(aws.head())

mr='sirivan_run4'
umhours=fio.model_outputs[mr]['filedates']
offset=timedelta(hours=10)
d0=umhours[0] + offset
dN=umhours[-1] + offset

AWS_plot_timeseries(aws,'FFDR/FFDI',d0=d0,dN=dN,color='r')
AWS_plot_timeseries(aws,'GFDR/GFDI',d0=d0,dN=dN,color='m')
ax=plt.gca()
plt.ylabel('Index')
ax2=ax.twinx()
AWS_plot_timeseries(aws,'T', d0=d0,dN=dN,color='k')
plt.ylabel("Degrees C")
plotting.add_legend(ax,colours=['r','m','k'],labels=['FFDI','GFDI','Temperature'])

cubes=fio.read_model_timeseries(mr,latlon=plotting._latlons_['uarbry'],dN=umhours[1])
#print(cubes)