# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 11:47:15 2019
    Examine AWS dataset(s)
    Currently only have Wagerup weather station data
    Compare model run outputs interpolated to lat/lon point where AWS sits
@author: jgreensl
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mpldates
from matplotlib import patches
import pandas
import iris
from copy import deepcopy
from datetime import datetime, timedelta

from utilities import fio, plotting, utils, constants

_sn_='AWS'

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
            ax.tick_params(axis='x',labelbottom=False)
        else:
            #plt.xlabel('time (AWST)')
            plt.xlabel(df.index.name)
        plt.title(title,fontsize=15)
        if units is not None:
            plt.ylabel(units[title],fontsize=12)
    return f, axes

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

def AIFS_Summary(mr='sirivan_run4',
                 d0=None,dN=None, d0_model=None, dN_model=None,
                 laptop=False):
    """
    show AIFS dataset sites, and compare T, RH, FFDI, and winds against 
    colocated model data
    ARGS:
        model_d0: datetime for beginning of model input (defualt full model output)
        model_dN: datetime for end of model input (defualt full model output)
        d0: datetime window left edge for plotting (default model time -2 hours)
        dN: datetime right edge (default model time + 2 hours)
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
    offset=timedelta(hours=fio.run_info[mr]['UTC_offset'])
                
    ## Subset to datetime using these limits
    if d0_model is None:
        d0_model=fio.run_info[mr]['filedates'][0]
    if dN_model is None:
        dN_model=fio.run_info[mr]['filedates'][-1]
    
    d0= d0_model+offset-timedelta(hours=3)
    dN= dN_model+offset+timedelta(hours=3)
                
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
        print("DEBUG: plotting ", sitei, site, latlon)
        color=sitecolors[sitei]
        # plot_date arguments for model timeseries
        pdargs = {'color':color,
                  'linestyle':ls_model,
                  'marker':''}
        # aws_timeseries arguments for site timeseries
        awsargs= {'d0':d0,
                  'dN':dN,
                  'color':color,
                  'linestyle':ls,
                  'linewidth':lw
                  }
        
        df = fio.AIFS_read_path(path)
        #site_wd=(df.loc[d0:dN])['Wd']
        #site_u,site_v = utils.uv_from_wind_degrees(site_wd.to_numpy())
        #site_lt = df.loc[d0:dN].index.to_numpy()
        #site_ws=(df.loc[d0:dN])['Ws km/h'].to_numpy() / 3.6 # km/h -> m/s
        
        
        cubes=None
        WESN=fio.run_info[mr]['WESN']
        # only have subset of data on laptop!
        if laptop:
            WESN=fio.run_info[mr]['WESN_laptop']
        latlon_model=np.copy(latlon)
        # if latlon is close to model region, compare nearby
        if latlon[1] < WESN[0]:
            if latlon[1] > WESN[0]-0.7:
                latlon_model=[latlon[0],WESN[0]+0.075]
            else:
                latlon_model=None
        elif latlon[1] > WESN[1]:
            if latlon[1] < WESN[1]+0.7:
                latlon_model=[latlon[0],WESN[1]-0.075]
            else:
                latlon_model=None
        if latlon[0] > WESN[3]:
            if latlon[0] < WESN[3]+0.6:
                if latlon_model is not None: 
                    latlon_model=[WESN[3]-0.075, latlon_model[1]]
                else:
                    latlon_model=[WESN[3]-0.075, latlon[1]]
            else:
                latlon_model=None
            
        if latlon_model is not None:
            print("DEBUG: WESN:",WESN)
            print("DEBUG: latlon, latlon_model:",latlon, latlon_model)
            print("DEBUG:", mr, d0_model, dN_model)
            
            cubes=fio.read_model_timeseries(mr,
                                            latlon=latlon_model,
                                            d0=d0_model,
                                            dN=dN_model,
                                            wind_10m=True)
            #print("DEBUG: === (timeseires) cubes ===")
            #print(cubes)
            # get model T, FDI, 10m Winds
            model_T,model_RH = cubes.extract(['air_temperature','relative_humidity'])
            model_ws,model_u,model_v = cubes.extract(['s_10m','u_10m','v_10m'])
            ctimes = utils.dates_from_iris(model_T) + offset
            model_T0 = model_T[:,0].data - 273.15 # Celcius
            model_RH0 = model_RH[:,0].data*100 # %
            model_ws0 = model_ws.data # m/s
            model_u0 = model_u.data # m/s
            model_v0 = model_v.data # m/s
            if np.any(model_ws0 < 0):
                print("ERROR: Windspeed < 0")
                print(model_ws0)
                print(model_u0)
                print(model_v0)
                print(np.sqrt(model_u0**2+model_v0**2))
                assert False, "Shouldn't have negative wind magnitude"
            Drought = np.repeat(10,len(model_T0))
            model_FFDI = utils.FFDI(Drought,model_RH0,model_T0,model_ws0*3.6)
            df_dti = pandas.DatetimeIndex(ctimes)
            model_df = pandas.DataFrame({'date':df_dti,'T':model_T0,
                                         'RH':model_RH0,'FFDI':model_FFDI,
                                         's':model_ws0})
            # model_df = pandas.DataFrame(ctimes,columns=['date'])
            model_df = model_df.set_index('date')
            print("DEBUG: model_df",model_df)
            print("DEBUG: site df:",df)
            
        
        ## Timeseries
        # Add time series for temperature
        plt.sca(ax_T)
        #print("DEBUG: df")
        #print(df)
        AWS_plot_timeseries(df,'T', **awsargs)
        if cubes is not None:
            #plt.plot_date(model_time_axis, model_T0, **pdargs)
            AWS_plot_timeseries(model_df, 'T', **pdargs)
        plt.ylabel('T [Celcius]')
        
        
        plt.sca(ax_RH)
        AWS_plot_timeseries(df,'RH', **awsargs)
        if cubes is not None:
            #plt.plot_date(model_time_axis,model_RH0, **pdargs)
            AWS_plot_timeseries(model_df, 'RH', **pdargs)
        plt.ylabel('RH [%]')
        
        # Add time series for indices
        plt.sca(ax_I)
        AWS_plot_timeseries(df,'FFDR/FFDI', **awsargs)
        if cubes is not None:
            # model FFDI
            #plt.plot_date(model_time_axis,model_FFDI, **pdargs)
            AWS_plot_timeseries(model_df, 'FFDI', **pdargs)
        #AWS_plot_timeseries(df,'GFDR/GFDI',d0=d0,dN=dN,color='m')
        plt.ylabel('FFDI')
        
        # add ts for winds
        plt.sca(ax_W)
        AWS_plot_timeseries(df,'Ws m/s', **awsargs)
        if cubes is not None:
            #plt.plot_date(model_time_axis, model_ws0, **pdargs)
            AWS_plot_timeseries(model_df,'s', **pdargs)
        plt.ylabel('Winds [m/s]')
        
        ## MAP ADDITIONS
        # Add point to map in axis 1
        plt.sca(ax1)
        plotting.map_add_nice_text(ax1,latlons=[latlon],texts=[site],markercolors=sitecolors[sitei])
        # add normal points to map
        plotting.map_add_locations_extent('sirivans',nice=True)
        if (cubes is not None) and (not np.all(np.array(latlon)==np.array(latlon_model))):
            plotting.map_add_nice_text(ax1,latlons=[latlon_model],texts=[''],markercolors=sitecolors[sitei], markers=['X'])
        
        plt.savefig("temp1.png")
    
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
    plt.xlabel("local time")
    #plt.xlim([d0,dN])
    fio.save_fig(mr,_sn_,"sirivan_AWS",plt)

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
                           groupbystr=None,
                           UTC=True):
    """
    read AWS dataframe, add model columns to dataframe
        optionally resample intervals (mean) using groupbystr
    
    arguments
    ---------
        groupbystr: string, make all columns binned ('30Min' by default)
        
    returns:
        df, subplots
            df: dataframe
            subplots: dictionary of metadata for plotting 
    """
    lat,lon = constants.latlons[AWS]
    extent = [lon-.02, lon+.02, lat-.02, lat+.02] # just grab real close to latlon
    heights=[1,10,30] # look at model heights 1m, 10m and 30m
    
    ## Read AWS: eventually handle multiple weather stations
    #if AWS == 'wagerup':
    subplots = deepcopy(__AWS_PLOTS__)
    df_aws, aws_attrs = fio.read_AWS_wagerup(UTC=UTC)

    umhours = fio.run_info[model_run]['filedates']
    
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
            # don't bother with a couple of metrics
            if (model_cube.name() == 'air_pressure') and (heights[zi]==30):
                continue
            elif (model_cube.name() == 's') and (heights[zi] == 1):
                continue
            
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
        heights = [1,10,30]
        hcolors = ['crimson','fuchsia','blueviolet']
        hmarkers = ['v','.','^']
        if mname=='s':
            heights = heights[1:]
            hcolors = hcolors[1:]
            hmarkers = hmarkers[1:]
        elif mname=='air_pressure':
            heights = heights[:2]
            hcolors = hcolors[:2]
            hmarkers = hmarkers[:2]
        subplots[title]['dfnames'].extend([mname+'%d'%d for d in heights])
        # add marker,color,linestyle too
        subplots[title]['colors'].extend(hcolors)
        subplots[title]['linestyles'].extend(['None']*len(hcolors))
        subplots[title]['markers'].extend(hmarkers)
    
    df = df_both
    if groupbystr is not None:
        df = df_both.resample('30Min')
    
    return df, subplots
    
def compare_site_to_model(AWS='wagerup',
                          model_run='waroona_run1',
                          showrange=None, # subset the time series to this range
                          UTC=True,
                          ):
    """
    AWS='wagerup', model_run='waroona_run1',dtoffset=timedelta(hours=8)):
    
    compare site to model run
    """
    
    df, subplots = combine_site_and_model(AWS=AWS, model_run=model_run, UTC=UTC)
    df_subset = df
    if showrange is not None:
        rangestr = [showrange[0].strftime("%Y-%m-%d %H"),
                     showrange[1].strftime("%Y-%m-%d %H")]
        df_subset = df.loc[rangestr[0]:rangestr[1]]
    
    ## now call method which plots the dataframe directed by my dictionary of subplot names
    df_time_series(df_subset, subplots=subplots)
    fio.save_fig(model_run=model_run, script_name=_sn_, 
                 pname='%s_%s_all%s'%(AWS,model_run,['','_UTC'][UTC]), plt=plt)
    
    ## aggregate and compare model to aws
    df_agg = df_subset.resample('30Min').mean()
    df_time_series(df_agg, subplots=subplots)
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
    fio.save_fig(model_run=model_run, script_name=_sn_, 
                 pname='%s_%s_correlation%s'%(AWS,model_run,['','_UTC'][UTC]), plt=plt)
    
    return df_agg, subplots
    
def plot_model_vs_data(df, colours=['maroon','k'], labels=['model','AWS']):
    """
    Show timeseries of model data superimposed on some other data
    
    INPUTS:
        df: dataframe containing: 
        dates: list of datetimes to be x axis
        colours: colours for model, meas
        labels: labels for model, meas
    """
    
    print("TODO")
    
    #df,subplots = AWS_series.combine_site_and_model()

    ### Site vs measurements 
    ## Show extent and weather station locations
    ## Read site and model using dataframe to aggregate on time:
    ## extract wd and ws arrays for comparison
    ## subplot for each site (argument list)

def HC06D_summary(mr='sirivan_run1'):
    """
    """
    # read site data
    data = fio.HC06D_read()
    #sites=pandas.concat(data.values(),axis=0)
    sitecolors=['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6',]
    ## Subset to datetime using these limits
    d0=fio.run_info[mr]['filedates'][0]
    dN=fio.run_info[mr]['filedates'][-1]
                
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
    
    fio.save_fig(mr,_sn_,"sirivan_AWS",plt)

if __name__=='__main__':
    
    ## Examine AIFS time series vs model weather/ffdi
    si_runs=['sirivan_run7_hr','sirivan_run5_hr','sirivan_run6_hr']
    
    #df = fio.AIFS_read_path(_AWS_['moree_airport']['path_AIFS'])
    #print(df)
    #df['T'].plot()
    
    d0=None
    dN=datetime(2017,2,11,22)#None
    if True:
        for run in si_runs:
            ## Test with dN so it doesn't take forever to read data
            #AIFS_Summary(run,d0=d0,dN_model=dN,laptop=True)
            AIFS_Summary(run)
    
    
    ## Summary of wagerup waroona stuff
    if False: 
        d0,dN = datetime(2016,1,4,10), datetime(2016,1,7,10)
        summary_wagerup(d0,dN,UTC=False)
        
        for mr in ['waroona_run2','waroona_run1','waroona_old']:
            for UTC in [False]:
                compare_site_to_model(AWS='wagerup',
                                      model_run=mr,
                                      showrange=[datetime(2016,1,5,6),datetime(2016,1,7,6)],
                                      UTC=UTC)


