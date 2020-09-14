#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 13:52:09 2019

  READING AND WRITING NETCDF AND GRIB(?) FILES

@author: jesse
"""

###
## IMPORTS
###

import iris
from iris.experimental.equalise_cubes import equalise_attributes
import numpy as np
import timeit # for timing stuff
import warnings
from datetime import datetime, timedelta
import pandas # for csv reading (AWS)

from glob import glob
import os

# This script is only run from parent folder, so this relative path should work
from utilities import utils

###
## GLOBALS
###
__VERBOSE__=True

## Info that is the same for most runs
#
sim_info={
    'sirivan':{
        'filedates':np.array([datetime(2017,2,11,21) + timedelta(hours=x) for x in range(24)]),
        'topog':'umnsaa_2017021121_slv.nc',
        'WESN':[148.9936, 150.6036, -32.8064, -31.1964],
        'UTC_offset':11, # utc + 11 in summer, utc+10 for other seasons
        },
    'waroona':{
        'filedates':np.array([datetime(2016,1,5,15) + timedelta(hours=x) for x in list(range(24))+list(range(36,60))]),
        'topog':'umnsaa_2016010515_slv.nc',
        'WESN':[115.2936,116.9036, -33.7064,-32.0964 ],
        'UTC_offset':8, # utc + 8, no summer daylight savings
        },
    'NYE':{
        'topog':'umnsaa_2019123009_slv.nc',
        'filedates':np.array([datetime(2019,12,30,9) + timedelta(hours=x) for x in range(24)]),
        'UTC_offset':11, # utc + 11 in summer
        }
    }


## Where are model outputs located, and run specific info
run_info = {
    ## run with new fire speed/fuel parameters changed again
    'sirivan_run7_hr':{
        'dir':'data/sirivan/run7/0p1/',
        'WESN':[149.3616,150.3975,-32.3024,-31.6985],
        'origdir':'/g/data/en0/hxy548/ACCESS-fire/sirivan/2020-08-19/20170211T2100Z/0p1/atmos/',
        },
    'sirivan_run7':{
        'dir':'data/sirivan/run7/0p3/',
        'WESN':[149.4412,150.2179,-32.3024,-31.6985],
        },
    ## run with new fire speed/fuel parameters turned up a bit again
    'sirivan_run6_hr':{
        'dir':'data/sirivan/run6/0p1/',
        'WESN':[149.4412,150.2179,-32.3024,-31.6985],
        'origdir':'/g/data/en0/hxy548/ACCESS-fire/sirivan/2020-08-19/20170211T2100Z/0p1/atmos/',
        },
    'sirivan_run6':{
        'dir':'data/sirivan/run6/0p3/',
        'WESN':[149.4412,150.2179,-32.3024,-31.6985],
        },
    ## run with new fire speed/fuel parameters turned up
    'sirivan_run5_hr':{
        'dir':'data/sirivan/run5/0p1/',
        'WESN':[149.2816,150.3175, -32.5184,-31.4825 ],
        },
    'sirivan_run5':{
        'dir':'data/sirivan/run5/0p3/',
        'WESN':[148.9936, 150.6036,-32.8064,-31.1964],
        'origdir':'/g/data/en0/hxy548/ACCESS-fire/sirivan/2020-08-14/20170211T2100Z/0p3/atmos/',
        },
    ## run with netcdf output (aug 4 2020)
    'sirivan_run4':{
        'dir':'data/sirivan/run4/0p3/',
        },
    ## New high res gadi run with fuelmap applied (April 2020)
    'sirivan_run3_hr':{
        'dir':'data/sirivan/run3_hr/',
        'fire_affix':'CSIRO_gadi_fuelmap',
        'origdir':'/scratch/en0/hxy548/cylc-run/au-aa860/share/cycle/20170211T2100Z/sirivan/0p1/ukv_os38/um/',
        'origfiredir':'/g/data/en0/hxy548/fire_vars/sirivan/0p1/'
        },
    ## New sirivan run (on GADI) by harvey in Feb 2020
    'sirivan_run2':{
        'dir':'data/sirivan/run2/',
        'fire_affix':'CSIRO_gadi',
        'origdir':'/scratch/en0/hxy548/cylc-run/au-aa860/share/cycle/20170211T2100Z/sirivan/0p3/ukv_os38/um/',
        'origfiredir':'/g/data/en0/hxy548/fire_vars/sirivan/0p3/'
        },
    ## high res version of sirivan run2 (100m x 100m)
    'sirivan_run2_hr':{
        'dir':'data/sirivan/run2_hr/',
        'fire_affix':'CSIRO_gadi',
        'run':'Run 26 Feb 2020',
        'origdir':'/scratch/en0/hxy548/cylc-run/au-aa860/share/cycle/20170211T2100Z/sirivan/0p1/ukv_os38/um/',
        'origfiredir':'/g/data/en0/hxy548/fire_vars/sirivan/0p1/'},
    ## New waroona run (on GADI) by harvey in Feb 2020
    ## Day 2 run on July 17th, 12 hour gap exists after the end of day 1,
    ## heat flux ramps up over 30 minutes, initial fuel load zero within fire ignition area
    'waroona_run3':{
        'dir':'data/waroona/run3/',
        'fire_affix':'CSIRO_new_gadi',
        'run':'Run 6 Feb 2020',
        'origdir':'/scratch/en0/hxy548/cylc-run/au-aa799/share/cycle/20160105T1500Z/waroona/0p3/ukv_os38/um/',
        'origdir2':'/scratch/en0/hxy548/cylc-run/au-aa876/share/cycle/20160107T0300Z/waroona/0p3/ukv_os38/um/',
        'origfiredir':'/g/data/en0/hxy548/fire_vars/waroona/0p3/'}, 
    ## New waroona run (on GADI) by harvey in July 2020
    ## same as waroona_run3 but fire coupling turned off
    'waroona_run3uc':{
        'dir':'data/waroona/run3uc/',
        'fire_affix':'CSIRO_gadi_ncp',
        'run':'Run 10 July 2020',
        'origdir':'/scratch/en0/hxy548/cylc-run/au-aa799/share/cycle/20160105T1500Z/waroona/0p3/ukv_os38/um/',
        'origdir2':'/scratch/en0/hxy548/cylc-run/au-aa876/share/cycle/20160107T0300Z/waroona/0p3/ukv_os38/um/',
        'origfiredir':'/g/data/en0/hxy548/fire_vars/waroona/0p3/'}, 
    ## Day 2 run in July 2020, Day 2 with no ramp for heat flux at fire ignition
    'waroona_run3_day2_orig':{
        'dir':'data/waroona/run3_day2_orig/',
        'filedates':np.array([datetime(2016,1,7,3) + timedelta(hours=x) for x in list(range(24))]),
        'fire_affix':'CSIRO_new_gadi',
        'run':'Run 6 Feb 2020',
        'origdir':'/scratch/en0/hxy548/cylc-run/au-aa876/share/cycle/20160107T0300Z/waroona/0p3/ukv_os38/um/',
        'origfiredir':'/g/data/en0/hxy548/fire_vars/waroona/0p3/'}, 
    ## Day 2 run in July, early start for day2 output with fire spreading too fast
    'waroona_run3_day2_early':{
        'dir':'data/waroona/run3_day2_early/',
        'filedates':np.array([datetime(2016,1,6,15) + timedelta(hours=x) for x in list(range(24))]),
        'fire_affix':'CSIRO_gadi',
        'run':'Run 13 July 2020',
        'origdir':'/scratch/en0/hxy548/cylc-run/au-aa876/share/cycle/20160106T1500Z/waroona/0p3/ukv_os38/um/',
        },
    ## Latest run from may 2020, 1p0 resolution run
    'waroona_run3_hr':{
        'dir':'data/waroona/run3_hr/',
        'run':'Run may 2020',
        'origdir':'/scratch/en0/hxy548/cylc-run/au-aa799/share/cycle/20160106T1500Z/waroona/1p0/ukv_os38/um/',
        },
    ## Copy of waroona_run2 suite, with fire coupling turned off
    'waroona_run2uc':{
        'dir':'data/waroona/run2uc/',
        'fireaffix':'fire/firefront.CSIRO_new_ncp.20160105T1500Z.nc',
        'path_fireflux':'fire/sensible_heat.CSIRO_new_ncp.20160105T1500Z.nc',
        'path_firespeed':'fire/fire_speed.CSIRO_new_ncp.20160105T1500Z.nc',
        'run':'Run ~12 December 2019',
        'origdir':'/short/en0/jwg574/cylc-run/au-aa799/share/cycle/20160105T1500Z/waroona/0p3/ukv_os38/um/',
        'origfiredir':'/raijin/short/en0/jwg574/cylc-run/au-aa799/work/20160105T1500Z/',
        },
    ## Run 1 was the first one I looked at, with east west rolls and no pyrocb
    ## Attemp to recreate 'old' run weather and pyrocb, updated fire to use mean met (rather than instantaneous)
    'waroona_run2':{
        'dir':'data/waroona/run2/',
        'fire_affix':'CSIRO_24h_new',
        'run':'Run ~10 December 2019',
        'origdir':'/short/en0/hxy548/cylc-run/au-aa799/share/cycle/20160105T1500Z/waroona/0p3/ukv_os38/um/',
        'origfiredir':'/short/en0/hxy548/tmp/waroona/0p3',
        },
    ## Run 1 was the first one I looked at, with east west rolls and no pyrocb
    'waroona_run1':{
        'dir':'data/waroona/run1/',
        'run':'Run in august 2019',
        'fire_affix':'CSIRO_24h',
        'origdir':'/short/en0/hxy548/cylc-run/au-aa799/share/cycle/20160105T1500Z/waroona/0p3/ukv_os38/um/',
        'origfiredir':'/short/en0/hxy548/tmp/waroona/0p3/',
        },
    ## Old run had pyrocb but also lots of high clouds and hooked F160 bases
    'waroona_old':{
        'dir':'data/waroona/old/',
        'filedates':np.array([datetime(2016,1,5,15) + timedelta(hours=x) for x in range(18)]),
        'run':'Run in August 2018',
        'firefront':'fire/firefront.01.nc', # under the data directory?
        'origdir':'/g/data1a/en0/mxp548/access-fire/waroona/run3/accessdata/',
        'origfiredir':'/short/en0/mxp548/cylc-run/au-aa714/work/20160105T1500Z/waroona_0p4_ukv_os38_um_fcst_000/',
        },
    ## Sirivan run (june 2019) at 100mx100m resolution
    'sirivan_run1_hr':{
        'dir':'data/sirivan/run1_hr/',
        'run':'Run in June 2019?',
        'fire_affix':'CSIRO_MinT',
        'origdir':'/short/en0/hxy548/cylc-run/au-aa860/share/cycle/20170211T2100Z/sirivan/0p1/ukv_os38/um',
        'origfiredir':'/short/en0/hxy548/cylc-run/au-aa860/work/20170211T2100Z/sirivan_0p1_ukv_os38_um_fcst_000/',
        },
    ## Sirivan run around June (2019)
    'sirivan_run1':{
        'dir':'data/sirivan/run1/',
        'topog':'umnsaa_pa2017021121.nc',
        'run':'Run in June 2019?',
        'fire_affix':'CSIRO_MinT',
        'origdir':'/short/en0/hxy548/cylc-run/au-aa860/share/cycle/20170211T2100Z/sirivan/0p3/ukv_os38/um',
        'origfiredir':'/short/en0/hxy548/cylc-run/au-aa860/work/20170211T2100Z/sirivan_0p3_ukv_os38_um_fcst_000/',
        },
    ## Old Old run has bad lat/lons...
    ## Note oldold run doesn't have QC at stage 5 resolution
    'waroona_oldold':{
        'dir':'data/waroona_oldold/',
        'topog':'stage5_sfc_orog.nc',
        'filedates':np.array([datetime(2016,1,5,15)]),
        'run':'Run in october 2016',
        'origdir':'/g/data1/en0/rjf548/fires/waroona.2016010615.vanj'
        },
    ## Nested run for NYE atmosphere around fire in SE NSW
    ## Uses newer model version: 11.6
    'NYE_run1':{
        'dir':'data/NYE_run1/',
        'run':'Run in June 2020',
        'origdir':'/home/574/jwg574/cylc-run/u-bq575/share/cycle/20191230T0900Z/SEForest/0p3/',
        },
    }
# populate run_info with sim_info
for run in run_info.keys():
    loc=run.split('_')[0]
    for k,v in sim_info[loc].items():
        # if sim location info not in run_info[run] dict, add it
        if k not in run_info[run].keys():
            run_info[run][k] = v

###
## METHODS
###

def _constraints_from_extent_(extent, constraints=None, tol = 0.0001):
    """
        Return iris constraints based on [WESN] lon,lat limits (with a tolerance)
        this only looks at the cell centres, so data with or without cell bounds compares equally
        additional constraints can be ampersanded to the lat lon constraints
    """
    West,East,South,North = extent
    # NEED TO COMPARE TO CELL MIDPOINT SO THAT BEHAVIOUR IS CONSTANT WHETHER OR NOT THERE ARE BOUNDS
    constr_lons = iris.Constraint(longitude = lambda cell: West-tol <= cell.point <= East+tol)
    constr_lats = iris.Constraint(latitude = lambda cell: South-tol <= cell.point <= North+tol)
    if constraints is not None:
        constraints = constraints & constr_lats & constr_lons
    else:
        constraints = constr_lats & constr_lons
    return constraints

def create_wind_profile_csv(model_run='waroona_run1'):
    # waroona latlon = -32.84, 115.93
    cubes = read_model_run(model_run,fdtime=None,
                           extent=[-32.83,115.92,-32.85,115.94], # just want tiny area at bottom of escarp
                           add_winds=True, add_z=True, add_topog=True,
                           add_RH=True)
    
    x,y,wdir = cubes.extract(['u','v','wind_direction'])
    z, topog = cubes.extract(['z_th','surface_altitude'])
    T, RH = cubes.extract(['air_temperature','relative_humidity'])
    
    # extract single lat lon
    cubes1 = iris.cube.CubeList([x,y,wdir,z,topog, T, RH])
    lat,lon = -32.84,115.93
    for i in range(len(cubes1)):
        cubes1[i] = cubes1[i].interpolate([('longitude',lon), ('latitude',lat)],
                                          iris.analysis.Linear())
    
    x1,y1,wdir1,z1,topog1,T1,RH1 = cubes1
    Tg = T1[:,0].data
    RHg = RH1[:,0].data
    wdg = wdir1[:,0].data
    wsg = np.sqrt(x1[:,0].data**2+y1[:,0].data**2)
    
    # height above ground level
    agl = z1.data - topog1.data
    agl_coord = iris.coords.AuxCoord(agl, var_name='height_agl', units='m')
    agl_coord.guess_bounds()
    
    # interpolate to our heights
    heights = np.arange(500,2001,500)
    cubes2 = iris.cube.CubeList([x1,y1,wdir1])
    ## interpolate to 500,1000,1500,2000m above ground level
    for i in range(len(cubes2)):
        # add agl coord
        cubes2[i].add_aux_coord(agl_coord,1) 
        cubes2[i] = cubes2[i].interpolate([('height_agl',heights)], iris.analysis.Linear())
    
    print(cubes2)
    
    ### Create csv with gridded horizontal wind speed and direction
    ## 
    ## Sample headers:
    headers=['Time',
             'Temperature',
             'Relative humidity',
             'Wind direction ground',
             'Wind direction 500m',
             'Wind direction 1000m',
             'Wind direction 1500m',
             'Wind direction 2000m',
             'Wind speed ground',
             'Wind speed 500m',
             'Wind speed 1000m',
             'Wind speed 1500m',
             'Wind speed 2000m']
    # time format (UTC)
    times = utils.dates_from_iris(cubes[0],remove_seconds=True)
    tstrings = [dt.strftime("%Y-%m-%dT%H:%M:%S")for dt in times]
    
    import csv
    [u2,v2,wd2c] = cubes2
    ws2 = np.sqrt(u2.data**2+v2.data**2)
    wd2 = wd2c.data
    with open('outputs/waroona_winds.csv',mode='w',newline='') as winds_file:
        fwriter=csv.writer(winds_file, delimiter=',')
        fwriter.writerow(headers)
        for i in range(len(times)):
            fwriter.writerow([tstrings[i],Tg[i]-273.15,RHg[i], 
                              wdg[i], wd2[i,0], wd2[i,1], wd2[i,2], wd2[i,3],
                              wsg[i], ws2[i,0], ws2[i,1], ws2[i,2], ws2[i,3]])
    print("INFO: SAVED FILE outputs/waroona_winds.csv")

def read_AWS_wagerup(UTC=True):
    '''
    Read AWS csv file into pandas dataframe and return the frame
    30m wind direction is modulated to [0,360] range (was [-360,360]??)
    units and descriptions are read too
    
    Data details:
    McCaw, L., Burrows, N., Beecham, B., and Rampant, P. (2016). Reconstruction of
    the spread and behaviour of the Waroona bushfire (Perth Hills 68) 6-7
    January 2016. Government of Western Australia Department of Parks and
    Wildlife.
    https://library.dbca.wa.gov.au/static/FullTextFiles/072096.pdf
    
    RETURNS: dataframe, dictionary

    '''
    wagerup_path = 'data/AWS/Wagerup.csv'
    ## Read header info
    wagerup_header = pandas.read_csv(wagerup_path, nrows=17)
    # just three columns wanted
    wagerup_header = wagerup_header.loc[:,['Sensor','Units','Description']]
    # dictionary from header info:
    wagerup_attrs = {}
    for colname, units, desc in zip(wagerup_header['Sensor'].values,
                                    wagerup_header['Units'].values,
                                    wagerup_header['Description'].values):
        # remove ampersand and semicolon, and remove whitespace
        units = units.replace('&deg;','deg').strip()
        # remove <sup></sup> tags
        units = units.replace('<sup>','').replace('</sup>','')
        wagerup_attrs[colname] = {'units':units, 'desc':desc}

    ## Read data
    wagerup = pandas.read_csv(wagerup_path,skiprows=19)

    # Read datetimes, WAST is Australian Western Standard Time: UTC+08:00
    awst_offset = timedelta(hours=8)
    dates0 = [datetime.strptime(dstr, "%d/%m/%Y %H:%M") for dstr in wagerup.loc[:,'WAST']]
    dates1 = [datetime.strptime(dstr, "%d/%m/%Y %H:%M")-awst_offset for dstr in wagerup.loc[:,'WAST']]
    # DATA APPEARS TO BE IN UTC, WAST IS MISLEADING
    wagerup_attrs['utc0'] = {'units':'%Y%m%d %H%M%S','desc':'first date in utc: %s'%dates0[0].strftime("%Y%m%d %H%M%S")}
    wagerup_attrs['local_time0'] = {'units':'%Y%m%d %H%M%S','desc':'first date in local time: %s'%dates1[0].strftime("%Y%m%d %H%M%S")}

    # convert column to datetime64 in the dataframe
    wagerup['WAST'] = pandas.to_datetime(wagerup['WAST'], dayfirst=True)
    
    # THIS DATASET APPEARS TO NOT BE IN WAST
    # DATETIME - 8 or +16 hours matches a midnight minima in solar radiation..
    # rename WAST to UTC and make it index
    # unless we want local time (then convert utc to local time)
    # Aus West Std Time = UTC+08
    # perhaps they did utc - 8 instead of +8??
    #   in this case, local time is inbuilt + 16, UTC is inbuild + 8, 
    # or else they did local time +8
    #   so local time is inbuilt - 8, UTC is inbuilt inbuilt-16
    # or maybe they just did UTC and called it WAST
    #   so local time is inbuilt + 8
    mistake1=False # this one appears reasonably against waroona_old
    mistake2=True # this one does not really match wind direction, but temperature is great
    mistake3=False # this one does not appear reasonable against waroona_old, model temperature is backwards
    if UTC:    
        wagerup = wagerup.rename(columns={'WAST':'UTC'})
        wagerup = wagerup.set_index('UTC')
        if mistake1:
            wagerup.index = wagerup.index + pandas.DateOffset(hours=8)
        elif mistake2:
            wagerup.index = wagerup.index + pandas.DateOffset(hours=-16)
    else:
        # set index to localtime
        wagerup = wagerup.set_index('WAST')
        # AWST = UTC + 8 hours # Aust west std time
        if mistake1:
            wagerup.index = wagerup.index + pandas.DateOffset(hours=16)
        elif mistake2:
            wagerup.index = wagerup.index + pandas.DateOffset(hours=-8)
        elif mistake3:
            wagerup.index = wagerup.index + pandas.DateOffset(hours=8)

    # add hour as column (maybe for collating later)
    #wagerup['Hour'] = wagerup.index.hour

    # select data for a single day or time window
    #day1 = wagerup.loc['2016-01-06']
    #some_hours = wagerup.loc['2016-01-06 23':'2016-01-07 12'] #includes end points

    # replace negative wind directions with their value plus 360 (for consistency)
    # np.sum(wagerup['Dta30'] < 0) # this happens 20 times from 1909 entries !?
    wagerup['Dta30'] = wagerup['Dta30'] % 360

    # Remove negative solar radiation
    wagerup.SR = wagerup.SR.astype(float)
    newSR = wagerup['SR'].values
    newSR[newSR<0] = np.NaN
    wagerup.SR.values[:] = newSR

    return wagerup, wagerup_attrs

def read_pft(model_run='waroona_run1',times=None, lats=None, lons=None):
    '''
    read pft from file, interpolate to lats,lons if desired
    returns pft, pftlats, pftlons
    '''
    # load the single cube in PFT.nc
    fpath = "data/PFT/"+model_run+".nc"
    pft0, = iris.load(fpath)
    
    # pft is [time, lats, lons]
    pft = pft0
    if lats is not None:
        pft = pft.interpolate([('latitude',lats)],
                                iris.analysis.Linear())
    if lons is not None:
        pft = pft.interpolate([('longitude',lons)],
                                iris.analysis.Linear())
    
    plats,plons = pft.coord('latitude').points, pft.coord('longitude').points
    ptimes = utils.dates_from_iris(pft, remove_seconds=True)
    pftd = pft.data
    if times is not None:
        tslice = np.array([utils.date_index(time, ptimes) for time in times])
        ptimes = ptimes[tslice]
        pftd = np.squeeze(pftd[tslice,:,:])
        
    return pftd, ptimes, plats, plons 

def read_nc_iris(fpath, constraints=None, keepvars=None, HSkip=None):
    '''
    Read netcdf file using iris, returning cubeslist
    actual data is not read until you call cube.data
    constraints can be applied, or variable names, or both
    '''

    print("INFO: Reading(iris) ",fpath)
    # First read all cubes in file using constraints
    if constraints is not None:
        cubes = iris.load(fpath, constraints)
    else:
        cubes = iris.load(fpath)

    # If just some variables are wanted, pull them out
    if keepvars is not None:
        cubes = cubes.extract(keepvars)

    # Maybe we want to cut down on the horizontal resolution
    if HSkip is not None:
        if not (HSkip == False):
            # For each cube, apply ::HSkip to lon/lat dimension
            small_cubes = iris.cube.CubeList()
            for cube in cubes:
                if cube.ndim == 2:
                    mini = cube[::HSkip,::HSkip]
                elif cube.ndim == 3:
                    mini = cube[...,::HSkip,::HSkip]
                elif cube.ndim == 4:
                    mini = cube[...,...,::HSkip,::HSkip]
                else:
                    print("ERROR: HSKIP is missing the cube: ")
                    print(cube)
                    assert False, "shouldn't miss any cubes"
                small_cubes.append(mini)
            cubes = small_cubes
    if len(cubes) < 1:
        print("ERROR: NO CUBES RETURNED IN fio.read_nc_iris()")
        print("     :",constraints)
    return cubes

def fire_paths(model_run):
    fdir = run_info[model_run]['dir']+'fire/'
    affix='20' # if there is no affix, look for name.YYYYMMDDTHHmmZ.nc
    if 'fire_affix' in run_info[model_run]:
        affix=run_info[model_run]['fire_affix']
    #print("DEBUG: looking for file:",fdir+'firefront.'+affix+'*')
    ffpaths=glob(fdir+'firefront.'+affix+'*')
    fluxpaths=glob(fdir+'sensible_heat.'+affix+'*')
    fspaths=glob(fdir+'fire_speed.'+affix+'*')
    v10paths=glob(fdir+'10m_vwind.'+affix+'*')
    u10paths=glob(fdir+'10m_uwind.'+affix+'*')
    ffpaths.sort()
    fluxpaths.sort()
    fspaths.sort()
    v10paths.sort()
    u10paths.sort()
    return [ffpaths, fluxpaths, fspaths, u10paths, v10paths]

def read_fire(model_run='waroona_run1',
              dtimes=None, constraints=None, extent=None,
              firefront=True,
              sensibleheat=False,
              firespeed=False,
              wind=False,
              day1=None,
              day2=False,
              HSkip=None):
    '''
    Read fire output cubes matching dtimes time dim
    output is transposed from [time,lon,lat] -> [time, lat, lon] to match model output
    ARGUMENTS:
    '''
    ## If no fire exists for model run, return None
    fdir = run_info[model_run]['dir']+'fire/'
    if not os.path.exists(fdir):
        print("ERROR: no such filepath:",fdir)
        # needs to be iterable to match cubelist return type (with wind counted twice) 
        return [None]*np.sum([firefront,sensibleheat,firespeed,wind,wind]) 
    
    ## if dtimes are passed, maybe we want some of day1 and day2,
    ## set day flag automatically using dtimes if available
    ## if day flags are already set then just move on
    if day1 is None:
        if dtimes is not None:
            day1,day2 = False,False
            if "waroona" in model_run and dtimes[-1] >= datetime(2016,1,6,15,1):
                day2 = True
            if "waroona" in model_run and dtimes[0] < datetime(2016,1,6,15,1):
                day1 = True
            if "sirivan" in model_run:
                day1 = True
                day2 = False
        else:
            day1=True
        #modelhours = run_info[model_run]['filedates']
        #print("DEBUG: fire_read recursion")
        #print("     : model hours 0, -1: ",modelhours[0], modelhours[-1])
        #print("     : model hours 25, 24: ",modelhours[25], modelhours[24])
        #print("     : dtimes 0, -1: ",dtimes[0],dtimes[-1])
        #print("     : day1, day2:", day1, day2)
    print("DEBUG: dayflags:",day1,day2)
    ## if reading both days, read one at a time and combine
    if day1 and day2:
        cubelist1=read_fire(
            model_run=model_run,
            dtimes=dtimes, 
            constraints=constraints,
            extent=extent,
            firefront=firefront,
            sensibleheat=sensibleheat,
            firespeed=firespeed,
            wind=wind,
            day1=True,
            day2=False,
            HSkip=HSkip,
            )
        cubelist2=read_fire(
            model_run=model_run,
            dtimes=dtimes, 
            constraints=constraints,
            extent=extent,
            firefront=firefront,
            sensibleheat=sensibleheat,
            firespeed=firespeed,
            wind=wind,
            day1=False,
            day2=True,
            HSkip=HSkip,
            )
        ret_cubes = iris.cube.CubeList()
        for [item1, item2] in zip(cubelist1, cubelist2):

            cubelist = iris.cube.CubeList([item1,item2])
            iris.util.unify_time_units(cubelist)
            iris.experimental.equalise_cubes.equalise_attributes(cubelist)
            items = cubelist.concatenate()
            if len(items) > 1:
                print(items)
                assert False, "Concatenate didn't work"
            ret_cubes.append(items[0])
        return ret_cubes

    ## otherwise read fire paths and return fire cubes
    #sensible_heat.[affix.]YYYYMMDDTHHmmZ.nc',
    #fire_speed, 10m_uwind, 10m_vwind
    fpathlists= fire_paths(model_run)
    
    dind=0
    if day2: 
        dind=1
    
    if extent is not None:
        constraints = _constraints_from_extent_(extent,constraints)

    # Build up cubelist based on which files you want to read
    cubelist = iris.cube.CubeList()
    flags = [firefront, sensibleheat, firespeed, wind, wind]
    units = [None, 'Watts/m2', 'm/s', 'm/s', 'm/s']
    for flag, paths, unit in zip(flags, fpathlists, units):
        if flag:
            if len(paths) < 1:
                print("ERROR: missing files")
                print("     :", fpathlists)
            cube, = read_nc_iris(paths[dind], constraints=constraints, HSkip=HSkip)
            if unit is not None:
                cube.units=unit
            cubelist.append(cube)

    # Fix old run time units
    if model_run == 'waroona_old':
        for cube in cubelist:
            cube.coord('time').units = 'seconds since 2016-01-05 15:00:00'
            cube.coord('latitude').units = 'degrees'
            cube.coord('longitude').units = 'degrees'
    
    if model_run == 'sirivan_run1':
        for cube in cubelist:
            cube.coord('time').units = 'seconds since 2017-02-11 21:00:00'
            # match coord system to um output
            cube.coord('latitude').coord_system=iris.coord_systems.GeogCS(6371229.0)
            cube.coord('longitude').coord_system=iris.coord_systems.GeogCS(6371229.0)

    # Subset by time if argument exists
    if dtimes is not None:
        for i in range(len(cubelist)):
            cubelist[i] = subset_time_iris(cubelist[i], dtimes)
    
    # finally put latitude before longitude
    for cube in cubelist:
        
        nx = len(cube.coord('longitude').points)
        if len(cube.shape) == 2:
            if cube.shape[0] == nx:
                cube.transpose()
        if len(cube.shape) == 3:
            if cube.shape[1] == nx:
                # t, lon, lat -> t, lat, lon
                cube.transpose([0, 2, 1])
        
    return cubelist

def read_model_run(model_version, fdtime=None, subdtimes=None, extent=None,
                   add_topog=True, add_z=False, add_winds=False, add_theta=False,
                   add_dewpoint=False, add_RH=False, HSkip=None, constraints=None):
    '''
    Read output from particular model run into cubelist, generally concatenates
    along the time dimension.

    INPUTS:
        model_version: string (see model_outputs keys)
        fdtime: datetime or datetimes, optional
            file output datetime[s] to be read, if None then read all files
        subdtimes: iterable of datetimes, optional
            after reading files, just keep these datetime slices
            None will return all datetimes
        extent: list, optional
            [West, East, South, North] lon,lon,lat,lat list used to spatially
            subset the data (if desired)
        add_<metric>: bool,
            add data cubes to returned cubelist
            WARNING: Some of these will use RAM
        HSkip: integer, optional
            slice horizontal dimension [::HSkip] to reduce resolution

    RETURNS:
        iris.cube.CubeList with standardised dimensions [time, level, latitude, longitude]
    
    0: air_pressure / (Pa)                 (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    1: air_pressure_at_sea_level / (Pa)    (time: 6; latitude: 14; longitude: 14)
    2: air_pressure_rho / (Pa)             (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    3: air_temperature / (K)               (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    4: cld_ice / (kg kg-1)                 (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    5: cld_water / (kg kg-1)               (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    6: specific_humidity / (kg kg-1)       (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    7: surface_air_pressure / (Pa)         (time: 6; latitude: 14; longitude: 14)
    8: surface_temperature / (K)           (time: 6; latitude: 14; longitude: 14)
    9: upward_air_velocity / (m s-1)       (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    10: x_wind / (m s-1)                    (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    11: y_wind / (m s-1)                    (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    # if cloud stuff is there
      : qc / (g kg-1)                       (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    ## add_topog
      : surface_altitude / (m)              (latitude: 14; longitude: 14)
    ## add_z
      : z_th / (m)                         (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    ## add_winds
      : u / (m s-1)                    (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
      : v / (m s-1)                    (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
      : s / (m s-1)                    (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
      : wind_direction / (degrees)          (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    ## add_dewpoint
      : vapour_pressure / (100 Pa)     (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
      : dewpoint_temperature / (K)     (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    ## add_RH
      : relative_humidity / (1)             (time: 6; model_level_number: 140; latitude: 14; longitude: 14)
    ## add_theta
      : potential_temperature / (K)                        (time ; model_level_number ; latitude ; longitude)
    '''

    ## make sure we have model run data
    assert model_version in run_info.keys(), "%s not yet supported by 'read_model_run'"%model_version

    ddir = run_info[model_version]['dir'] + 'atmos/'
    locname=model_version.split('_')[0]
    
    timelesscubes=[]
    allcubes=None

    ## No ftimes? set to all ftimes
    if fdtime is None:
        fdtime = run_info[model_version]['filedates']
            
    # make sure it's iterable
    if not hasattr(fdtime,'__iter__'):
        fdtime = [fdtime]
    fdtimes = np.array(fdtime)

    ## First read the basics, before combining along time dim

    ### SPECIFIC TO OLD RUN
    if model_version=='waroona_old':
        for dtime in fdtimes:
            cubelist = read_waroona_old(dtime)
            if allcubes is None:
                allcubes = cubelist
            else:
                allcubes.extend(cubelist)
            
    elif model_version=='sirivan_run1':
        for dtime in fdtimes:
            cubelist = read_sirivan_run1(dtime, HSkip=HSkip)
            if allcubes is None:
                allcubes = cubelist
            else:
                allcubes.extend(cubelist)
        
    elif '_run' in model_version:
        for dtime in fdtimes:
            slv,ro1,th1,th2 = read_standard_run(dtime, 
                                                mv=model_version, 
                                                HSkip=HSkip,
                                                constraints=constraints)

            # Remove specific humidity from slv, since we get the whole array from th1
            if len(slv.extract('specific_humidity')) > 0:
                slv.remove(slv.extract('specific_humidity')[0])

            if allcubes is None:
                allcubes=slv
            else:
                allcubes.extend(slv)
            # model output on rho levels and theta levels are slightly different
            # rename the rho levels air pressure (or else it's repeated)
            ro1.extract('air_pressure')[0].rename('air_pressure_rho')

            allcubes.extend(ro1)
            allcubes.extend(th1)
            allcubes.extend(th2)

        if add_z:
            # also read model height above ground level
            # This field is only in the first output file
            height_varname = 'height_above_reference_ellipsoid'
            if constraints is not None:
                height_varname = height_varname & constraints
            height_date = run_info[model_version]['filedates'][0]
            ro1_height_path = ddir+height_date.strftime('umnsaa_%Y%m%d%H_mdl_ro1.nc')
            zro, = read_nc_iris(ro1_height_path,
                                constraints = height_varname, 
                                HSkip=HSkip)
            th1_height_path = ddir+height_date.strftime('umnsaa_%Y%m%d%H_mdl_th1.nc')
            zth, = read_nc_iris(th1_height_path,
                                constraints = height_varname,
                                HSkip=HSkip)
            zro.rename('z_rho')
            zth.rename('z_th')
            if 1 in zth.shape:
                zth = iris.util.squeeze(zth)
            if 1 in zro.shape:
                zro = iris.util.squeeze(zro)
            timelesscubes.append(zth)
            timelesscubes.append(zro)

    elif model_version == "waroona_oldold":
        #0: surface_altitude / (m)              (latitude: 88; longitude: 88)
        #1: eastward_wind / (m s-1)             (t: 8; Hybrid height: 70; latitude: 88; longitude: 89)
        #2: northward_wind / (m s-1)            (t: 8; Hybrid height: 70; latitude: 88; longitude: 88)
        #3: upward_air_velocity / (m s-1)       (t: 8; Hybrid height: 70; latitude: 88; longitude: 88)
        #4: z_th / (m)                          (Hybrid height: 70; latitude: 88; longitude: 88)
        #5: z_rho / (m)                         (Hybrid height: 70; latitude: 88; longitude: 88)
        oldoldcubes = read_waroona_oldold()

        # rename height dim to match other runs
        for i in range(len(oldoldcubes)):
            if len(oldoldcubes[i].shape) > 2:
                nz = len(oldoldcubes[i].coord('Hybrid height').points)
                oldoldcubes[i].coord('Hybrid height').rename('level_height')
                # add model_level_number as aux_coord
                zdim = [0,1][len(oldoldcubes[i].shape)>3]
                mlncoord = iris.coords.Coord(np.arange(nz),
                                             standard_name='model_level_number')
                oldoldcubes[i].add_aux_coord(mlncoord, data_dims=zdim)
            

        # rename winds
        oldoldcubes.extract('eastward_wind')[0].rename('x_wind')
        oldoldcubes.extract('northward_wind')[0].rename('y_wind')

        allcubes = oldoldcubes.extract(['x_wind','y_wind','upward_air_velocity'])
        # rename time dim to 'time'
        for i in range(len(allcubes)):
            allcubes[i].coord('t').rename('time')

        # topog and z levels added automatically for oldold run (easier)
        timelesscubes.extend(oldoldcubes.extract(['surface_altitude','z_th','z_rho']))

    ### GENERIC TO ALL MODEL OUTPUT (some exceptions)

    ## Concatenate along time dimension
    ## First need to unify time dimension:
    iris.util.unify_time_units(allcubes)
    ## Also need to equalise the attributes list
    # I think this just deletes attributes which are not the same between matching cubes..
    equalise_attributes(allcubes)
    ## Join along the time dimension
    allcubes = allcubes.concatenate()
    ## subset to our subdtimes
    if subdtimes is not None:
        cubedates=utils.dates_from_iris(allcubes[0])
        tslice = np.array([utils.nearest_date_index(dt, cubedates, allowed_seconds=120) for dt in subdtimes])
        for i in range(len(allcubes)):
            allcubes[i] = allcubes[i][tslice]
    
    # surface stuff may have no time dim...
    for varname in ['surface_air_pressure','surface_temperature']:
        saps = allcubes.extract(varname)
        if len(saps) > 1:
            sap = saps.merge_cube()
            for cube in saps:
                allcubes.remove(cube)
            tube = allcubes[0]
            time = tube.coord('time')
            # add time steps for easyness
            sap0 = sap
            if tube.shape[0] != sap.shape[0]:
                sap0 = sap.interpolate([('time',time.points)],
                                        iris.analysis.Linear())
            # add single cube with time dim
            allcubes.append(sap0)
            
    ## NOW add any timeless cubes
    allcubes.extend(timelesscubes)

    # topography
    if add_topog:
        # oldold is special case, topog read in read_waroona_oldold()
        if model_version != "waroona_oldold":
            topog, = read_nc_iris(ddir+run_info[model_version]['topog'],
                                  constraints = 'surface_altitude',
                                  HSkip=HSkip)
            if 1 in topog.shape:
                topog = iris.util.squeeze(topog) # don't want time dim in topog
            allcubes.append(topog) 

    ## Subset spatially
    # based on extent
    if extent is not None:
        constraints = _constraints_from_extent_(extent)
        for i in range(len(allcubes)):
            allcubes[i] = allcubes[i].extract(constraints)

    ## extras
    # Rename some stuff in model versions 4 and above (new netcdf)
    # also remove dupes here
    #print("DEBUG:",allcubes)
    for dupecubes in ["mass_fraction_of_cloud_ice_in_air","upward_air_velocity","mass_fraction_of_cloud_liquid_water_in_air","air_temperature","air_pressure",]:
        dupes = allcubes.extract(dupecubes)
        if len(dupes) > 1:
            remove_duplicate_cubes(dupes,allcubes,ndims_expected=4)
    #if int(model_version[-1]) > 3:
    #    allcubes.extract('vertical_wnd')[0].rename('upward_air_velocity')
    
    # Add cloud parameter at g/kg scale
    if model_version not in ["sirivan_run1",]:
        # Rename ice and water cubes
        #if int(model_version[-1]) < 4:
        allcubes.extract('mass_fraction_of_cloud_liquid_water_in_air')[0].rename('cld_water')
        allcubes.extract('mass_fraction_of_cloud_ice_in_air')[0].rename('cld_ice')
        
        water_and_ice = allcubes.extract(['cld_water',
                                          'cld_ice'])
        if len(water_and_ice) == 2:
            water,ice=water_and_ice
            qc = (water+ice) * 1000
            qc.units = 'g kg-1'
            qc.rename('qc')
            allcubes.append(qc)

    if add_z:
        # This is down here so that it only happens AFTER subsetting
        if (model_version == 'waroona_old') or (model_version == 'sirivan_run1'):
            # add zth cube
            p, pmsl = allcubes.extract(['air_pressure','air_pressure_at_sea_level'])
            # take out time dimension
            p, pmsl = p[0], pmsl[0]
            nz,ny,nx = p.shape
            # repeat surface pressure along z axis
            reppmsl = np.repeat(pmsl.data[np.newaxis,:,:],nz, axis=0)
            zth = -(287*300/9.8)*np.log(p.data/reppmsl)
            iris.std_names.STD_NAMES['z_th'] = {'canonical_units': 'm'}
            zthcube=iris.cube.Cube(zth, standard_name='z_th',
                                   var_name="zth", units="m",
                                   dim_coords_and_dims=[(p.coord('model_level_number'),0),
                                                        (p.coord('latitude'),1),
                                                        (p.coord('longitude'),2)])
            allcubes.append(zthcube)

    if add_winds:
        # wind speeds need to be interpolated onto non-staggered latlons
        u1, v1 = allcubes.extract(['x_wind','y_wind'])
        
        ### DESTAGGER u and v using iris interpolate
        ### (this will trigger the delayed read)
        # u1: [t,z, lat, lon1]
        # v1: [t,z, lat1, lon]  # put these both onto [t,z,lat,lon]
        u = u1.interpolate([('longitude',v1.coord('longitude').points)],
                           iris.analysis.Linear())
        v = v1.interpolate([('latitude',u1.coord('latitude').points)],
                           iris.analysis.Linear())
        
        # add standard names for these altered variables:
        iris.std_names.STD_NAMES['u'] = {'canonical_units': 'm s-1'}
        iris.std_names.STD_NAMES['v'] = {'canonical_units': 'm s-1'}
        u.standard_name='u'
        v.standard_name='v'
        # Get wind speed cube using hypotenuse of u,v
        s = utils.wind_speed_from_uv_cubes(u,v)
        s.units = 'm s-1'
        s.var_name='s' # s doesn't come from a var with a std name so can just use var_name
        
        # Get wind direction using arctan of y/x
        wind_dir = utils.wind_dir_from_uv(u.data,v.data)
        
        wdcube = iris.cube.Cube(wind_dir,
                                var_name='wind_direction',
                                units='degrees',
                                #long_name='degrees clockwise from due North',
                                dim_coords_and_dims=[(s.coord('time'),0),
                                                     (s.coord('model_level_number'),1),
                                                     (s.coord('latitude'),2),
                                                     (s.coord('longitude'),3)])
        wdcube.units = 'degrees'
        wdcube.var_name='wind_direction'
        
        # add cubes to list
        allcubes.append(u)
        allcubes.append(v)
        allcubes.append(s)
        allcubes.append(wdcube)

    if add_dewpoint:
        # Take pressure and relative humidity
        #print("DEBUG: add_dewpoint", allcubes)
        p,q = allcubes.extract(['air_pressure','specific_humidity'])
        p_orig_units = p.units
        q_orig_units = q.units
        p.convert_units('hPa')
        q.convert_units('kg kg-1')
        #print("DEBUG: units", p_orig_units, p.units, q_orig_units,q.units)
        # calculate vapour pressure:
        epsilon = 0.6220 # gas constant ratio for dry air to water vapour
        e = p*q/(epsilon+(1-epsilon)*q)
        e.rename('vapour_pressure')
        e.units = 'hPa'
        #e_correct=np.mean(e.data)
        p.convert_units(p_orig_units)
        q.convert_units(q_orig_units)
        #assert np.isclose(np.mean(e.data),e_correct), "Changing units back messes with vapour_presssure"

        allcubes.append(e)
        # calculate dewpoint from vapour pressure
        Td = 234.5 / ((17.67/np.log(e.data/6.112))-1) # in celcius
        Td = Td + 273.15 # celcius to kelvin

        # change Td to a Cube
        iris.std_names.STD_NAMES['dewpoint_temperature'] = {'canonical_units': 'K'}
        cubeTd = iris.cube.Cube(Td, standard_name="dewpoint_temperature",
                                   var_name="Td", units="K",
                                   dim_coords_and_dims=[(p.coord('time'),0),
                                                        (p.coord('model_level_number'),1),
                                                        (p.coord('latitude'),2),
                                                        (p.coord('longitude'),3)])

        allcubes.append(cubeTd)

    if add_theta:
        # Estimate potential temp
        p, Ta = allcubes.extract(['air_pressure','air_temperature'])
        theta = utils.potential_temperature(p.data,Ta.data)
        # create cube
        iris.std_names.STD_NAMES['potential_temperature'] = {'canonical_units': 'K'}
        cubetheta = iris.cube.Cube(theta, standard_name="potential_temperature",
                                   var_name="theta", units="K",
                                   dim_coords_and_dims=[(p.coord('time'),0),
                                                        (p.coord('model_level_number'),1),
                                                        (p.coord('latitude'),2),
                                                        (p.coord('longitude'),3)])
        allcubes.append(cubetheta)
    
    if add_RH:
        # estimate relative humidity
        q,T = allcubes.extract(['specific_humidity','air_temperature'])
        # compute RH from specific and T in kelvin
        orig_Tunits=T.units
        T.convert_units('K')
        RH = utils.relative_humidity_from_specific(q.data.data, T.data.data)
        # restore T units (just in case)
        T.convert_units(orig_Tunits)
        # turn RH into a cube and add to return list
        iris.std_names.STD_NAMES['relative_humidity'] = {'canonical_units': '1'}
        cubeRH = iris.cube.Cube(RH, standard_name="relative_humidity",
                                   var_name="RH", units="1",
                                   dim_coords_and_dims=[(q.coord('time'),0),
                                                        (q.coord('model_level_number'),1),
                                                        (q.coord('latitude'),2),
                                                        (q.coord('longitude'),3)])
        allcubes.append(cubeRH)
    
    return allcubes


def read_standard_run(dtime, constraints=None, extent=None, mv='waroona_run1', HSkip=None):
    '''
        Read converted model output files
        returns list of 4 iris cube lists, matching the model output files: slv, ro1, th1, th2
        ========
        0: specific_humidity / (1)             (time: 6; latitude: 576; longitude: 576)
        1: surface_air_pressure / (Pa)         (time: 6; latitude: 576; longitude: 576)
        2: air_pressure_at_sea_level / (Pa)    (time: 6; latitude: 576; longitude: 576)
        3: surface_temperature / (K)           (time: 6; latitude: 576; longitude: 576)
        ========
        0: air_pressure / (Pa)                 (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        1: x_wind / (m s-1)                    (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        2: y_wind / (m s-1)                    (time: 6; model_level_number: 140; latitude: 577; longitude: 576)
        ========
        0: air_pressure / (Pa)                 (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        1: air_temperature / (K)               (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        2: specific_humidity / (kg kg-1)       (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        3: upward_air_velocity / (m s-1)       (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        ========
        0: mass_fraction_of_cloud_ice_in_air / (kg kg-1) (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        1: mass_fraction_of_cloud_liquid_water_in_air / (kg kg-1) (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
    '''

    dstamp=dtime.strftime('%Y%m%d%H')
    ddir = run_info[mv]['dir']+'atmos/'

    # If we just want a particular extent, subset to that extent using constraints
    if extent is not None:
        constraints = _constraints_from_extent_(extent,constraints)

    # 3 file types for new waroona output
    ## slx - single level output
    ## mdl_ro1 - multi level, on ro dimension (winds)
    ## mdl_th1 - multi level, on theta dimension
    ## mdl_th2 - multi level, on theta dimension (since um output file size is limited)
    _standard_run_vars_ = {
        'slv':[
            'specific_humidity',  # kg/kg [t,lat,lon]
            'spec_hum', # kg/kg [t,lat,lon] # new runs
            'surface_air_pressure', # Pa
            'air_pressure_at_sea_level', # Pa [t,lat,lon]
            'surface_temperature', # [t,y,x]
            # x_wind and y_wind are here also...
            ],
        'mdl_ro1':[
            'air_pressure', # [t,z,lat,lon]
            'x_wind', # [t,z,y,x]
            'y_wind', # [t,z,y,x]
            #'height_above_reference_ellipsoid', # m [z,y,x]
            ],
        'mdl_th1':[
            'air_pressure', # Pa [t,z,y,x]
            'pressure', # Pa [tzyx] # new runs
            'air_temperature', # K [t,z,y,x]
            'air_temp', # K [tzyx] # this one is comes from new runs
            'specific_humidity', # kg/kg [tzyx]
            'upward_air_velocity', # m/s [tzyx]
            'vertical_wnd', # m/s [tzyx]
            ],
        'mdl_th2':[
            'mass_fraction_of_cloud_ice_in_air', # kg/kg [tzyx]
            'mass_fraction_of_cloud_liquid_water_in_air',
            'cld_water',
            'cld_ice',
            ],
        }
    cubelists = []
    for filetype,varnames in _standard_run_vars_.items():
        path = ddir+'umnsaa_%s_%s.nc'%(dstamp,filetype)
        cubelists.append(read_nc_iris(path,constraints=constraints,keepvars=varnames,HSkip=HSkip))

    return cubelists

def read_waroona_old(dtime, constraints=None, extent=None):
    '''
    0: mass_fraction_of_cloud_liquid_water_in_air / (kg kg-1) (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    1: mass_fraction_of_cloud_ice_in_air / (kg kg-1) (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    2: upward_air_velocity / (m s-1)       (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    3: air_pressure / (hPa)                (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    4: air_temperature / (K)               (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    5: x_wind / (m s-1)                    (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    6: y_wind / (m s-1)                    (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    7: air_pressure_at_sea_level / (Pa)    (time: 2; latitude: 88; longitude: 106)
    8: specific_humidity / (kg kg-1)       (time: 2; model_level_number: 140; latitude: 88; longitude: 106)

    '''
    dstamp=dtime.strftime('%Y%m%d%H')
    ddir = run_info['waroona_old']['dir']+'atmos/'
    # First read topography

    # If we just want a particular extent, subset to that extent using constraints
    if extent is not None:
        constraints = _constraints_from_extent_(extent,constraints)

    path = ddir+'umnsaa_pc%s.nc'%dstamp
    varnames = ['mass_fraction_of_cloud_liquid_water_in_air',
                'mass_fraction_of_cloud_ice_in_air',
                'upward_air_velocity', # m/s [t,z,lat,lon]
                'air_pressure', # Pa [t,z,lat,lon]
                'surface_air_pressure', # Pa [lat,lon] ?? no time dim?
                'air_temperature', # K [lat,lon]
                'surface_temperature', # K [lat,lon] ?? no time dim??
                'air_temperature_0', # K [t, z, lat, lon]
                'x_wind','y_wind', # m/s [t,z,lat,lon]
                'air_pressure_at_sea_level', # Pa [time, lat, lon]
                'specific_humidity', # kg/kg [lat,lon]
                'specific_humidity_0', # kg/kg [t,z,lat,lon] #???
                ]
    cubes = read_nc_iris(path,constraints=constraints,keepvars=varnames)

    # specific_humidity is the name of two or three variables, we just want the good one
    sh_cubes = cubes.extract('specific_humidity')
    remove_duplicate_cubes(sh_cubes,cubes)
    # same deal with air temp: seems to be issue with iris cubes package?
    ta_cubes = cubes.extract('air_temperature')
    remove_duplicate_cubes(ta_cubes,cubes)

    return cubes

def read_waroona_oldold(constraints=None, extent=None):
    '''
        read waroona_oldold model output
        0: surface_altitude / (m)              (latitude: 88; longitude: 88)
        1: eastward_wind / (m s-1)             (t: 8; Hybrid height: 70; latitude: 88; longitude: 89)
        2: northward_wind / (m s-1)            (t: 8; Hybrid height: 70; latitude: 88; longitude: 88)
        3: upward_air_velocity / (m s-1)       (t: 8; Hybrid height: 70; latitude: 88; longitude: 88)
        4: z_th / (m)                          (Hybrid height: 70; latitude: 88; longitude: 88)
        5: z_rho / (m)                         (Hybrid height: 70; latitude: 88; longitude: 88)
    '''
    ddir = run_info['waroona_oldold']['dir']
    xwind_path = ddir+'combined_alltimes_ml_xwind_stage5.nc'
    ywind_path = ddir+'combined_alltimes_ml_ywind_stage5.nc'
    zwind_path = ddir+'combined_alltimes_ml_zwind_stage5.nc'

    #113.916, 113.9208, 113.9256, ... 118.5096, 118.5144, 118.5192
    lons = np.linspace(113.916, 118.5192, 960, endpoint=True)
    lons1 = np.linspace(113.9184, 118.5216, 960, endpoint=True) # staggered
    # -35.73, -35.726, -35.722, -35.718, ... -30.942, -30.938, -30.934
    lats = np.linspace(-35.73, -30.934, 1200, endpoint=True)
    lats1 = np.linspace(-35.728, -30.936, 1199, endpoint=True) # staggered

    # set up dimension coordinates
    latdim = iris.coords.DimCoord(lats,'latitude')
    latdim1 = iris.coords.DimCoord(lats1,'latitude')
    londim = iris.coords.DimCoord(lons,'longitude')
    londim1 = iris.coords.DimCoord(lons1,'longitude')

    cubelist=iris.cube.CubeList()
    if os.path.isfile(ddir+'oldold_xwind_s5_subset.nc'):
        print("INFO: on local machine, reading subset data")
        xwind_path = ddir+'oldold_xwind_s5_subset.nc'
        ywind_path = ddir+'oldold_ywind_s5_subset.nc'
        zwind_path = ddir+'oldold_zwind_s5_subset.nc'

    if extent is not None:
        constraints = _constraints_from_extent_(extent,constraints)

    ## Read data
    with warnings.catch_warnings():
        # ignore UserWarning: Ignoring netCDF variable 'hybrid_ht' invalid units 'level'
        warnings.simplefilter('ignore')
        xwind1, = read_nc_iris(xwind_path)
        ywind1, = read_nc_iris(ywind_path)
        zwind, = read_nc_iris(zwind_path)
        ## Read topography
        topog1, = read_nc_iris(ddir+run_info['waroona_oldold']['topog'])
        ## Read heights of theta and rho levels
        zth1, = read_nc_iris(ddir+'stage5_ml_htheta.nc')
        zrho1, = read_nc_iris(ddir+'stage5_ml_hrho.nc')

    # length one time dim and levels dim removed
    topog = topog1[0,0]
    zth   = zth1[0]
    zrho   = zrho1[0]
    
    # update z_theta and z_rho names
    iris.std_names.STD_NAMES['z_th'] = {'canonical_units': 'm'}
    iris.std_names.STD_NAMES['z_rho'] = {'canonical_units': 'm'}
    zth.standard_name='z_th'
    zrho.standard_name='z_rho'

    # add latlon dims to data
    topog.add_dim_coord(latdim,0)
    topog.add_dim_coord(londim,1)

    xwind1.add_dim_coord(latdim,2)
    xwind1.add_dim_coord(londim1,3)

    ywind1.add_dim_coord(latdim1,2)
    ywind1.add_dim_coord(londim,3)

    zwind.add_dim_coord(latdim,2)
    zwind.add_dim_coord(londim,3)

    for cube in [zth, zrho]:
        cube.add_dim_coord(latdim,1)
        cube.add_dim_coord(londim,2)

    # add cube to list (subset if possible)
    for cube in [topog, xwind1, ywind1, zwind, zth, zrho]:
        if constraints is not None:
            cube=cube.extract(constraints)
        cubelist.append(cube)

    return cubelist

def read_sirivan_run1(dtime, constraints=None, extent=None, HSkip=None):
    '''
    0: mass_fraction_of_cloud_liquid_water_in_air / (kg kg-1) (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    1: mass_fraction_of_cloud_ice_in_air / (kg kg-1) (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    2: upward_air_velocity / (m s-1)       (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    3: air_pressure / (hPa)                (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    4: air_temperature / (K)               (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    5: x_wind / (m s-1)                    (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    6: y_wind / (m s-1)                    (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    7: air_pressure_at_sea_level / (Pa)    (time: 2; latitude: 88; longitude: 106)
    8: specific_humidity / (kg kg-1)       (time: 2; model_level_number: 140; latitude: 88; longitude: 106)

    '''
    dstamp=dtime.strftime('%Y%m%d%H')
    ddir = run_info['sirivan_run1']['dir']+'atmos/'
    # First read topography

    # If we just want a particular extent, subset to that extent using constraints
    if extent is not None:
        constraints = _constraints_from_extent_(extent,constraints)

    path = ddir+'umnsaa_pc%s.nc'%dstamp
    varnames = ['mass_fraction_of_cloud_liquid_water_in_air',
                'mass_fraction_of_cloud_ice_in_air',
                'upward_air_velocity', # m/s [t,z,lat,lon]
                'air_pressure', # Pa [t,z,lat,lon]
                'air_temperature', # K [lat,lon] ?
                'air_temperature_0', # K [t, z, lat, lon]
                'x_wind','y_wind', # m/s [t,z,lat,lon]
                'air_pressure_at_sea_level', # Pa [time, lat, lon]
                'specific_humidity', # kg/kg [lat,lon]
                'specific_humidity_0', # kg/kg [t,z,lat,lon] #???
                ]
    cubes = read_nc_iris(
        path,
        constraints=constraints,
        keepvars=varnames,
        HSkip=HSkip
        )

    # specific_humidity is the name of two or three variables, we just want the good one
    sh_cubes = cubes.extract('specific_humidity')
    remove_duplicate_cubes(sh_cubes,cubes,ndims_expected=4)
    # same deal with air temp: seems to be issue with iris cubes package?
    ta_cubes = cubes.extract('air_temperature')
    remove_duplicate_cubes(ta_cubes,cubes,ndims_expected=4)

    return cubes

def remove_duplicate_cubes(dupes, cubeslist, ndims_expected=3):
    flag = False
    if len(dupes) > 1:
        for dupe in dupes:
            if flag or (len(dupe.shape) != ndims_expected):
                cubeslist.remove(dupe)
            else:
                # keep first 3d instance of sh
                flag = True

def read_topog(model_version, extent=None, HSkip=None):
    '''
    Read topography cube
    '''

    ddir = run_info[model_version]['dir']+'atmos/'

    constraints='surface_altitude'
    if extent is not None:
        constraints = _constraints_from_extent_(extent,constraints)
    topog, = read_nc_iris(ddir + run_info[model_version]['topog'],
                          constraints = constraints, HSkip=HSkip)

    # don't want time dim in topog
    topog = iris.util.squeeze(topog) 
    return topog

def subset_time_iris(cube,dtimes,seccheck=121):
    '''
    take a cube with the time dimension and subset it to just have dtimes
    can handle iris seconds, minutes, or hours time dim formats
    assert times are available within seccheck seconds of dtimes
    '''
    tdim  = cube.coord('time')
    secmult = 1
    grain = str(tdim.units).split(' ')[0]
    unitformat = '%s since %%Y-%%m-%%d %%H:%%M:00'%grain
    if grain == 'minutes':
        secmult=60
    elif grain == 'hours':
        secmult=3600

    d0 = datetime.strptime(str(tdim.units),unitformat)
    # datetimes from ff
    dt = np.array([d0 + timedelta(seconds=secs*secmult) for secs in tdim.points])
    # for each datetime in dtimes argument, find closest index in dt
    tinds = []
    for dtime in dtimes:
        tinds.append(np.argmin(abs(dt-dtime)))
    tinds = np.array(tinds)

    # Check that fire times are within 2 minutes of desired dtimes
    #print("DEBUG: diffs")
    #print([(ffdt[tinds][i] - dtimes[i]).seconds < 121 for i in range(len(dtimes))])
    assert np.all([(dt[tinds][i] - dtimes[i]).total_seconds() < seccheck for i in range(len(dtimes))]), "fire times are > 2 minutes from requested dtimes"

    # subset cube to desired times
    return cube[tinds]

def make_folder(pname):
    folder = '/'.join(pname.split('/')[:-1]) + '/'
    if not os.path.exists(folder):
        print("INFO: Creating folder:",folder)
        os.makedirs(folder)

def save_fig_to_path(pname,plt, **savefigargs):
    '''
    Create dir if necessary
    Save figure
    example: save_fig('my/path/plot.png',plt)
    INPUTS:
        pname = path/to/plotname.png
        plt = matplotlib.pyplot instance

    '''
    # Defaults:
    if 'dpi' not in savefigargs:
        savefigargs['dpi']=150
    if 'transparent' not in savefigargs:
        savefigargs['transparent']=False
        
    make_folder(pname)
    print ("INFO: Saving figure:",pname)
    plt.savefig(pname, **savefigargs)
    plt.close()

def standard_fig_name(model_run, script_name, pname, 
                      subdir=None, ext='.png'):
    if isinstance(pname,datetime):
        dstamp=pname.strftime('%Y%m%d%H%M')
        pname='fig_%s%s'%(dstamp,ext)

    if subdir is not None:
        pname='%s/'%subdir + pname

    if len(pname.split('.')) == 1: # no extension
        pname=pname+ext
    
    path='figures/%s/%s/%s'%(model_run,script_name,pname)
    return path

def save_fig(model_run, script_name, pname, plt, subdir=None,
             ext='.png', **savefigargs):
    """
    create figurename as figures/model_run/script_name/[subdir/]fig_YYYYMMDDhhmm.png

    INPUTS:
        model_run,script_name : strings
        pname : can be datetime or string
            if datetime, pname is set to fig_YYYYMMDDhhmm
        plt : is instance of matplotlib.pyplot
        subdir : string - optional extra subdir
        ext : optional, only used if pname has no extension
    """
    path = standard_fig_name(model_run, script_name, pname, subdir, ext)

    save_fig_to_path(path,plt, **savefigargs)


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

## Read AIFS based csv files
def AIFS_read_path(path):
    """
    Return pandas dataframe based on path input
    titles:
        Station	DF	Cur %	Grass t/ha	EST HH:MM dd/mm/yy	
        T	Td	RH	Wd	Ws km/h	Wg km/h	FFDR/FFDI	GFDR/GFDI
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

def read_model_timeseries(model_run,latlon,
                          radial_avg_degrees=None,
                          d0=None,dN=None,
                          wind_10m=True,
                          ):
    """
    ARGUMENTS:
        radial_avg_degrees: instead of interpolating to a lat/lon, average horizontally up to this many degrees away
            currently this is a square not a circle averaged
        d0,dN: start and end date (optional)
        wind_10m flag tells method to add 10m winds from fire output files
    """
    lat,lon = latlon
    extent = [lon-.02, lon+.02, lat-.02, lat+.02] # just grab real close to latlon
    if radial_avg_degrees is not None:
        #assert radial_avg_degrees>0.01, "radial average too small!"
        extent = [lon-radial_avg_degrees, 
                  lon+radial_avg_degrees,
                  lat-radial_avg_degrees, 
                  lat+radial_avg_degrees]
    
    umhours = run_info[model_run]['filedates']
    
    # limit hours desired to d0-dN
    if d0 is not None:
        di=min(utils.date_index(d0,umhours)-1, 0)
        umhours=umhours[di:]
    if dN is not None:
        di=utils.date_index(dN,umhours)
        umhours=umhours[:di+1]
    ## Read Model output:
    cubes = read_model_run(model_run, fdtime=umhours, extent=extent, 
                                 add_topog=True, add_winds=True, 
                                 add_RH=True, add_z=True, 
                                 add_dewpoint=True, add_theta=True)
    
    ctimes = utils.dates_from_iris(cubes.extract('upward_air_velocity')[0])
    
    if wind_10m:
        u10,v10=read_fire(model_run, extent=extent, dtimes=ctimes,
                          firefront=False, wind=True)
        u10.rename('u_10m')
        v10.rename('v_10m')
        s10=utils.wind_speed_from_uv_cubes(u10,v10)
        s10.units='m s-1'
        s10.var_name='s_10m'
        cubes.append(u10)
        cubes.append(v10)
        cubes.append(s10)
    
    # remove unwanted times
    if d0 is not None:
        di=utils.date_index(d0,ctimes)
        for ci,cube in enumerate(cubes):
           if 'time' in cube.coords():
               cubes[ci] = cube[di:]
        ctimes=ctimes[di:]
    if dN is not None:
        di=utils.date_index(dN,ctimes)
        for ci,cube in enumerate(cubes):
           if 'time' in cube.coords():
               cubes[ci] = cube[:(di+1)]
        ctimes=ctimes[:di]
    # Interpolate or average the horizontal component
    for ci,cube in enumerate(cubes):
        #print("DEBUG: interp step:",ci,cube.name(), cube.shape)
        if 'latitude' in [coord.name() for coord in cube.coords()]:
            if radial_avg_degrees:
                # dimensions to collapse will be 'latitude' and 'longitude'
                cubes[ci] = cube.collapsed(['latitude','longitude'], iris.analysis.MEAN)
            else:
                cubes[ci] = utils.profile_interpolation(cube,latlon,)
        #print("     :",ci,cubes[ci].name(),cubes[ci].shape)
    return cubes
