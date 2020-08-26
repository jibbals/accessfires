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

_sn_='localtests'


# Dict to hold all the AWS metadata
# Normal resolution run bounds:
#lat_min = -32.8064 ;
#lat_max = -31.1964 ;
#lon_min = 148.9936 ;
#lon_max = 150.6036 
# High res run 5 bounds:
#lat_min = -32.5184 ;
#lat_max = -31.4825 ;
#lon_min = 149.2816 ;
#lon_max = 150.3175 ;
# high res run after run6
#lat_min = -32.3024 ;
#lat_max = -31.6985 ;
#lon_min = 149.4412 ;
#lon_max = 150.2179 ;


_AWS_={
    "moree_airport":{
        "path_AIFS":"data/AWS/MoreeAirport.csv",
        "latlon":[-29.4946, 149.8505], # degrees
        "latlon_model":None , # no equivalent nearby spot within model bounds
        "latlon_model_hr":None , 
        "altitude":214, # metres ASL
        },
    "coonabarabran_airport":{
        "path_AIFS":"data/AWS/CoonabarabranAirport.csv",
        "latlon":[-31.3330,149.2699],
        "latlon_model":[-31.3330,149.2699],
        "latlon_model_hr":[-31.51,149.2699],
        "altitude":645,
        },
    "dubbo_airport":{
        "path_AIFS":"data/AWS/DubboAirport.csv",
        "latlon":[-32.2189,148.5696],
        "latlon_model":[-32.22,149.1], 
        "latlon_model_hr":[-32.22,149.3], 
        "altitude":285,
        },
    "murrurundi_gap":{
        "path_AIFS":"data/AWS/MurrurundiGap.csv",
        "latlon":[-31.74,150.79],
        "latlon_model":[-31.74,150.5],
        "latlon_model_hr":[-31.74,150.3],
        "altitude":729.4,
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
    # add column for Ws in m/s
    # km/h * 1000/3600
    data['Ws m/s'] = data['Ws km/h'] / 3.6
    
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

def AIFS_Summary(mr='sirivan_run5',d0=None,dN=None):
    """
    show AIFS dataset sites, and compare T, RH, FFDI, and winds against 
    colocated model data
    """
    sites=[]
    for site in _AWS_.keys():
        if 'path_AIFS' in _AWS_[site]:
            sites.append(site)
    sitecolors=['#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6',]
    # linestyles for model/data
    ls_model='--'
    ls='-'
    lw=3
    offset=timedelta(hours=10)
                
    ## Subset to datetime using these limits
    h0_model=fio.model_outputs[mr]['filedates'][0]
    hN_model=fio.model_outputs[mr]['filedates'][-1]
    d0= h0_model+offset-timedelta(hours=2)
    dN= hN_model+offset+timedelta(hours=2)
                
    # create figure    
    fig=plt.figure(figsize=[13,20])
    # first subplot is map
    _,ax1 = plotting.map_tiff_qgis(fname='sirivan_AIFS.tiff',fig=fig,subplot_row_col_n=[5,1,1])
    
    # Timeseries: T, indices, Wind
    ax_T = plt.subplot(5,1,2)
    ax_RH= plt.subplot(5,1,3,sharex=ax_T)
    ax_I = plt.subplot(5,1,4,sharex=ax_T)
    ax_W = plt.subplot(5,1,5,sharex=ax_T)
    
    # for each site, get latlon, color, plot data
    for sitei,site in enumerate(sites):
        path = _AWS_[site]['path_AIFS']
        latlon = _AWS_[site]['latlon']
        
        color=sitecolors[sitei]
        
        pdargs = {'color':color,
                  'linestyle':ls_model,
                  'marker':''}
        awsargs= {'d0':d0,
                  'dN':dN,
                  'color':color,
                  'linestyle':ls,
                  'linewidth':lw
                  }
        quiverargs={'color':color, 
                    'alpha':0.5,
                    'pivot':'mid', # put arrow on centre of ws line
                    #'scale':0.05, # scale down the auto size??
                    'headwidth':3,
                    'headlength':2,
                    'headaxislength':2,
                    'width':.004
                    }
        
        df = AIFS_read_path(path)
        site_wd=(df.loc[d0:dN])['Wd']
        site_u,site_v = utils.uv_from_wind_degrees(site_wd.to_numpy())
        site_lt = df.loc[d0:dN].index.to_numpy()
        site_ws=(df.loc[d0:dN])['Ws km/h'].to_numpy() / 3.6 # km/h -> m/s
        
        
        cubes=None
        latlon_model = _AWS_[site]['latlon_model']
        if '_hr' in mr:
            latlon_model = _AWS_[site]['latlon_model_hr']
        if latlon_model is not None:
            #print("WARNING: using uarbry instead of site latlon until on NCI")
            #cubes=fio.read_model_timeseries(mr, 
            #                                latlon=plotting._latlons_['uarbry'], 
            #                                dN=h0_model+timedelta(hours=2),
            #                                wind_10m=True)
            #print("DEBUG: ",cubes)
            print("INFO: reading", latlon_model)
            cubes=fio.read_model_timeseries(mr,
                                            latlon=latlon_model,
                                            dN=hN_model,
                                            wind_10m=True)
            #print(cubes)
            # get model T, FDI, 10m Winds
            model_T,model_RH = cubes.extract(['air_temperature','relative_humidity'])
            model_ws,model_u,model_v = cubes.extract(['s_10m','u_10m','v_10m'])
            ctimes = utils.dates_from_iris(model_T) + offset
            model_T0 = model_T[:,0].data.data - 273.15
            model_RH0 = model_RH[:,0].data*100 # %
            model_ws0 = model_ws.data # m/s
            model_u0 = model_u.data # m/s
            model_v0 = model_v.data # m/s
            Drought = np.repeat(10,len(model_T0))
            model_FFDI = utils.FFDI(Drought,model_RH0,model_T0,model_ws0*3.6)
        
        ## Timeseries
        # Add time series for temperature
        plt.sca(ax_T)
        AWS_plot_timeseries(df,'T', **awsargs)
        if cubes is not None:
            plt.plot_date(ctimes, model_T0, **pdargs)
        #plt.title('T')
        plt.ylabel('T [Celcius]')
        
        plt.sca(ax_RH)
        AWS_plot_timeseries(df,'RH', **awsargs)
        if cubes is not None:
            plt.plot_date(ctimes,model_RH0, **pdargs)
        plt.ylabel('RH [%]')
        
        # Add time series for indices
        plt.sca(ax_I)
        AWS_plot_timeseries(df,'FFDR/FFDI', **awsargs)
        if cubes is not None:
            # model FFDI
            plt.plot_date(ctimes,model_FFDI, **pdargs)
        #AWS_plot_timeseries(df,'GFDR/GFDI',d0=d0,dN=dN,color='m')
        plt.ylabel('FFDI')
        
        # add ts for winds
        plt.sca(ax_W)
        AWS_plot_timeseries(df,'Ws m/s', **awsargs)
        # now add quiver for sites
        n_arrows=12
        qskip = max(int(np.floor(len(site_u)/n_arrows)),1)
        als = 0.6 # arrow length scale
        plt.quiver(site_lt[::qskip], site_ws[::qskip], 
                   site_u[::qskip]*als, site_v[::qskip]*als, 
                   **quiverargs,
                   )
        
        if cubes is not None:
            plt.plot_date(ctimes, model_ws0, **pdargs)
        plt.ylabel('Winds [m/s]')
    
        if cubes is not None:
            # normalize windspeed for unit length quivers
            wdnorm = np.sqrt(model_u0**2 + model_v0**2)
            # dont show quiver at every single point
            qskip = max(int(np.floor(len(wdnorm)/n_arrows)),1)
            # Add quiver
            plt.quiver(ctimes[::qskip], model_ws0[::qskip], 
                       model_u0[::qskip]/wdnorm[::qskip]*als, 
                       model_v0[::qskip]/wdnorm[::qskip]*als, 
                       **quiverargs,
                       )
        
        ## MAP ADDITIONS
        # Add point to map in axis 1
        plt.sca(ax1)
        plotting.map_add_nice_text(ax1,latlons=[latlon],texts=[site],markercolors=sitecolors[sitei])
        # add normal points to map
        plotting.map_add_locations_extent('sirivans',nice=True)
        if (cubes is not None) and (not np.all(np.array(latlon)==np.array(latlon_model))):
            plotting.map_add_nice_text(ax1,latlons=[latlon_model],texts=[''],markercolors=sitecolors[sitei], markers=['X'])
    
    # legend for sites
    plotting.add_legend(ax_T,
                        colours=sitecolors[:len(sites)],
                        labels=sites,
                        )
    # legend for markers
    plotting.add_legend(ax_I,
                        colours=['k']*2,
                        labels=['data','model'],
                        styles=['-','--'],
                        )
    plt.sca(ax_W)
    plt.xlabel("local time (UTC+10)")
    fio.save_fig(mr,_sn_,"sirivan_AWS",plt)

#AIFS_Summary('sirivan_run5_hr')

### READ MOREE
#sname='moree_airport'
#path = _AWS_[sname]['path_AIFS']
#latlon = _AWS_[sname]['latlon']
#
#aws=AIFS_read_path(path)
#print(aws.head())
#site_u,site_v = utils.uv_from_wind_degrees(aws['Wd'].to_numpy())
#site_lt = aws.index.to_numpy()
#ws=aws['Ws km/h'].to_numpy()
#plt.quiver(site_lt, ws , site_u, site_v)
#
#
mr='sirivan_run5_hr'
umhours=fio.model_outputs[mr]['filedates'][:2]
cubes=fio.read_model_run(mr,extent=[150.0,150.2,-31.99,-31.9],fdtime=umhours)
print(cubes)

T=cubes.extract('air_pressure')[0]
h0,h1 = utils.height_from_iris(T,bounds=True)[0]
print("DEBUG: h0, h1:",h0,h1)
tmean=T[0,25:48,:,:].collapsed('model_level_number',iris.analysis.MEAN)
h0,h1 = utils.height_from_iris(tmean,bounds=True)[0]
print("DEBUG: meaned h0, h1:",h0,h1)

