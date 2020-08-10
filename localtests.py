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

HC06D_numlatlon={
        #'COONABARABRAN':['064008',-31.2786,149.2786], # DAILY
        'COONABARABRAN AIRPORT':['064017',-31.3330,149.2699],
        #'DUNEDOO':['064009',-32.0165,149.3956], # DAILY!
        #'GULGONG':['062013',-32.3634,149.5329], # DAILY!
        'MERRIWA':['061287',-32.19,150.17],
        'MUDGEE':['062101',-32.5628,149.6149],
        'NULLO MOUNTAIN':['062100',-32.7244,150.2290],
        #'WELLINGTON':['065034',-32.5635,148.9503], # DAILY!
        }

def HC06D_read_path(path):
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
    
    ## Read data
    print("INFO: reading ",path)
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
    # Remove white space from station name:
    data['Station Name'] = data['Station Name'].str.strip()
    # only keep the columns we want
    data = data.loc[:,HC06D_dict.values()]
    
    # Check all quality flags are "N"ormal (I'm assuming that's what N is for)
    qkeys = []
    for key in HC06D_dict.values():
        if 'Q_' in key:
            # indices of bad data
            bads = data[key] != "N"
            # if any are bad, reset to NaN
            if np.sum(bads) > 0:
                # Q_<index_name>, want to fix the index_name column
                data.loc[bads,key[2:]] = np.NaN
                print("Warning:", np.sum(bads), "entries removed from", key[2:], "at %s"%data["Station Name"][0] )
            qkeys.append(key)
    
    # finally drop quality flag column
    data = data.drop(columns=qkeys)
    
    # make sure data columns are numeric (not string)
    for key in data.keys():
        if key not in ['Station Name','LocalTime']:
            data[key] = pandas.to_numeric(data[key], errors='coerce')
    
    # set local time to be the index
    data = data.set_index('LocalTime')
    
    return data
    
def HC06D_read(sites=[], extent=None, concat=False):
    """
    READ AVAILABLE HC06D sites
    ARGUMENTS:
        subset to sites list (optional)
        subset to extent [left,right,bot,top] in degrees (optional)
        concat combines the sites into one dataframe
    """
    
    # path template
    pbase = 'data/AWS/HC06D_Data_%s_48628999806477.txt' 
    
    ret_data={} # dataframes to be returned
    if len(sites) > 0:
        # subset by site names
        for site in sites:
            path=pbase%(HC06D_numlatlon[site][0])
            ret_data[site]=HC06D_read_path(path)
    elif extent is not None:
        # subset by extent
        for site in HC06D_numlatlon.keys():
            num,lat,lon = HC06D_numlatlon[site]
            if (lat>extent[2]) and (lat<extent[3]) and (lon>extent[0]) and (lon<extent[1]):
                path=pbase%(num)
                ret_data[site]=HC06D_read_path(path)
    else:
        # grab all
        for site in HC06D_numlatlon.keys():
            path=pbase%(HC06D_numlatlon[site][0])
            ret_data[site]=HC06D_read_path(path)
    
    if concat:
        return pandas.concat(ret_data.values(),axis=0)
    return ret_data
    #dates = [datetime.strptime(dstr, "%d/%m/%Y %H:%M") for dstr in data.loc[:,dind_str]]
    
    
    
#### Show site locations/temperatures/winds etc...

# read site data
data = HC06D_read()
#sites=pandas.concat(data.values(),axis=0)
sitecolors=['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6',]
## Subset to datetime using these limits
d0=fio.model_outputs['sirivan_run1']['filedates'][0]
dN=fio.model_outputs['sirivan_run1']['filedates'][-1]
            
# create figure    
fig=plt.figure(figsize=[13,16])
# first subplot is map
_,ax1 = plotting.map_tiff_qgis(fname='sirivans.tiff',fig=fig,subplot_row_col_n=[3,1,1])

for itemi,itemname in enumerate(['Ta','ws']):
    axi=plt.subplot(3,1,itemi+2)
    latlons=[]
    for sitei,site in enumerate(data.keys()):
        #print("INFO: plotting",key)
        df=data[site].loc[d0.strftime("%Y-%m-%d %H:%M"):dN.strftime("%Y-%m-%d %H:%M")]
        df[itemname].plot(color=sitecolors[sitei])
        latlons.append([data[site]['Latitude'][0],data[site]['Longitude'][0]])
    plt.title(itemname)
    
# Add point to map in axis 1
plotting.map_add_nice_text(ax1,latlons=latlons,texts=data.keys(),markercolors=sitecolors[:len(data.keys())])
# add normal points to map
plt.sca(ax1)
plotting.map_add_locations_extent('sirivans',nice=True)

plotting.add_legend(axi,colours=sitecolors[:len(data.keys())],labels=data.keys())


#for i,sitename in enumerate(HC06D_numlatlon.keys()):
#    # grab site, subset by datetime
#    df = data[sitename].loc[d0.strftime("%Y-%m-%d %H:%M"):dN.strftime("%Y-%m-%d %H:%M")]
#    dates=df.index.values
#    plt.plot_date(dates,df['Ta'],color=sitecolors[i],label=sitename)
#plt.legend()
#

