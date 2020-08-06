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

def read_HC06D_path(path):
    """
    Return Pandas table with Columns based on following 
        original names:table names
    {   'Day/Month/Year Hour24:Minutes in DD/MM/YYYY HH24:MI format in Local Time':'LocalTime',
        'Air temperature in Degrees C':'Ta',
        'Quality of air temperature':'Q_Ta',
        'Dew point temperature in Degrees C':'Td',
        'Quality of dew point temperature':'Q_Td',
        'Wet bulb temperature in Degrees C':'Twb',
        'Quality of wet bulb temperature':'Q_Twb',
        'Relative humidity in percentage %':'RH',
        'Quality of relative humidity':'Q_RH',
        'Wind speed measured in m/s':'ws',
        'Quality of wind speed':'Q_ws',
        'Wind direction measured in degrees':'wd',
        'Quality of wind direction':'Q_wd',
        'Mean sea level pressure in hPa':'Pmsl',
        'Quality of mean sea level pressure':'Q_Pmsl',
        'Station level pressure in hPa':'Pstn',
        'Quality of station level pressure':'Q_Pstn',
        'QNH pressure in hPa':'Pqnh',
        'Quality of QNH pressure':'Q_Pqnh',
        'Vapour pressure in hPa':'Pv',
        'Quality of vapour pressure':'Q_Pv',
        'Saturated vapour pressure in hPa':'Psv',
        'Quality of saturated vapour pressure':'Q_Psv',
        'Total cloud amount in eighths':'Cloud',
        'Quality of total cloud amount':'Q_Cloud',
    }
    """
    path='data/AWS/HC06D_Data_061287_48628999806477.txt'
    ## Read data
    data = pandas.read_csv(path) # skiprows=1 removes titles

    # read datetimes
    dind_str='Day/Month/Year Hour24:Minutes in DD/MM/YYYY HH24:MI format in Local Time'
    # 01/01/2017 00:00
    # convert column to datetime64 in the dataframe
    data[dind_str] = pandas.to_datetime(data[dind_str], dayfirst=True)
    
    HC06D_dict = {
        'Station Name':'Station Name',
        'Locality':'Locality',
        'Latitude':'Latitude',
        'Longitude':'Longitude',
        'Height above MSL':'Height above MSL',
        'Day/Month/Year Hour24:Minutes in DD/MM/YYYY HH24:MI format in Local Time':'LocalTime',
        'Air temperature in Degrees C':'Ta',
        'Quality of air temperature':'Q_Ta',
        'Dew point temperature in Degrees C':'Td',
        'Quality of dew point temperature':'Q_Td',
        'Wet bulb temperature in Degrees C':'Twb',
        'Quality of wet bulb temperature':'Q_Twb',
        'Relative humidity in percentage %':'RH',
        'Quality of relative humidity':'Q_RH',
        'Wind speed measured in m/s':'ws',
        'Quality of wind speed':'Q_ws',
        'Wind direction measured in degrees':'wd',
        'Quality of wind direction':'Q_wd',
        'Mean sea level pressure in hPa':'Pmsl',
        'Quality of mean sea level pressure':'Q_Pmsl',
        'Station level pressure in hPa':'Pstn',
        'Quality of station level pressure':'Q_Pstn',
        'QNH pressure in hPa':'Pqnh',
        'Quality of QNH pressure':'Q_Pqnh',
        'Vapour pressure in hPa':'Pv',
        'Quality of vapour pressure':'Q_Pv',
        'Saturated vapour pressure in hPa':'Psv',
        'Quality of saturated vapour pressure':'Q_Psv',
        'Total cloud amount in eighths':'Cloud',
        'Quality of total cloud amount':'Q_Cloud',
    }
    
    # rename using dict
    data = data.rename(columns=HC06D_dict)
    # set local time to be the index
    data = data.set_index('LocalTime')
    return data.loc[:,HC06D_dict.values()]
    
def read_HC06D(sites=[], extent=None):
    """
    READ AVAILABLE HC06D sites
    ARGUMENTS:
        subset to sites list (optional)
        subset to extent [left,right,bot,top] in degrees (optional)
    """
    site_numlatlon={
        'COONABARABRAN':['064008',-31.2786,149.2786],
        'COONABARABRAN AIRPORT':['064017',-31.3330,149.2699],
        'DUNEDOO':['064009',-32.0165,149.3956],
        'GULGONG':['062013',-32.3634,149.5329],
        'MERRIWA':['061287',-32.19,150.17],
        'MUDGEE':['062101',-32.5628,149.6149],
        'NULLO MOUNTAIN':['062100',-32.7244,150.2290],
        'WELLINGTON':['065034',-32.5635,148.9503],
        }
    #site_numlatlon['NULLO MOUNTAIN']=site_numlatlon['NULLO'] # alternate name
    
    
    ret_data={} # dataframes to be returned
    if len(sites) > 0:
        # subset by site names
        for site in sites:
            path='HC06D_Data_%s_48628999806477.txt'%(site_numlatlon[site][0])
            ret_data[site]=read_HC06D_path(path)
    elif extent is not None:
        # subset by extent
        for site in site_numlatlon.keys():
            num,lat,lon = site_numlatlon[site]
            if (lat>extent[2]) and (lat<extent[3]) and (lon>extent[0]) and (lon<extent[1]):
                path='HC06D_Data_%s_48628999806477.txt'%(num)
                ret_data[site]=read_HC06D_path(path)
    else:
        # grab all
        for site in site_numlatlon.keys():
            path='HC06D_Data_%s_48628999806477.txt'%(site_numlatlon[site][0])
            ret_data[site]=read_HC06D_path(path)
    return ret_data
    #dates = [datetime.strptime(dstr, "%d/%m/%Y %H:%M") for dstr in data.loc[:,dind_str]]