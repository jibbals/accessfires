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
import numpy as np
import timeit # for timing stuff
from datetime import datetime, timedelta

from glob import glob
import os


### 
## GLOBALS
###

model_outputs = {
        # Attemp to recreate 'old' run weather and pyrocb
        'waroona_run2':{
                'path':'data/waroona_run2/',
                'topog':'',
                'run':'Run in September 2019',
                'origdir':'/short/en0/hxy548/tmp/waroona/0p3/',
                'origfiredir':''},
        # Run 1 was the first one I looked at, with east west rolls and no pyrocb
        'waroona_run1':{
                'path':'data/waroona_run1/',
                'topog':'umnsaa_2016010515_slv.nc',
                'run':'Run in august 2019',
                'origdir':'/short/en0/hxy548/cylc-run/au-aa799/share/cycle/20160105T1500Z/waroona/0p3/ukv_os38/um/',
                'origfiredir':'/short/en0/hxy548/tmp/waroona/0p3/'},
        # Old run had pyrocb but also lots of high clouds and hooked F160 bases
        'waroona_old':{
                'path':'data/waroona_old/',
                'topog':'',
                'run':'Run in August 2018',
                'origdir':'/g/data1a/en0/mxp548/access-fire/waroona/run3/accessdata'},
        # Old Old run has bad lat/lons...
        'waroona_oldold':{
                'path':'data/waroona_oldold/',
                'topog':'stage5_sfc_orog.nc',
                'run':'Run in october 2016',
                'origdir':'/g/data1/en0/rjf548/fires/waroona.2016010615.vanj'},
                 }


def save_fig(pname,plt):
    '''
    Create dir if necessary
    Save figure
    example: save_fig('my/path/plot.png',plt)
    INPUTS:
        pname = path/to/plotname.png
        plt = matplotlib.pyplot instance 
    
    '''
    folder = '/'.join(pname.split('/')[:-1]) + '/'
    if not os.path.exists(folder):
        print("INFO: Creating folder:",folder)
        os.makedirs(folder)
    print ("INFO: Saving figure:",pname)
    plt.savefig(pname)
    plt.close()

#from .context import utils
# Couldn't import utils, just copying this function for now
def potential_temperature(p,T):
    '''
    calculate theta from pressure and air temperature
    # Potential temperature based on https://en.wikipedia.org/wiki/Potential_temperature
    # with gas constant R = 287.05 and specific heat capacity c_p = 1004.64
    '''
    nt,nz,ny,nx = p.shape
    Ta  = T[:,0:1,:,:] # [t,1,lat,lon] at surface
    repTa = np.repeat(Ta[:,:,:,:], nz, axis=1) # repeat Ta along z axis
    assert np.all(repTa[:,0,:,:] - repTa[:,1,:,:] == 0), "Repeated z dim is not the same"
    return repTa*(1e5/p)**(287.05/1004.64)

###
## GLOBALS
###
__VERBOSE__=True

## Sir ivan pc fire files
_topog_sirivan_ = 'data/sirivan/umnsaa_pa2017021121.nc'
_topog_waroona_ = 'data/waroona/umnsaa_pa2016010515.nc'
_files_sirivan_ = sorted(glob('data/sirivan/umnsaa_pc*.nc'))


# 3 file types for new waroona output
## slx - single level output
## mdl_ro1 - multi level, on ro dimension (winds)
## mdl_th1 - multi level, on theta dimension
## mdl_th2 - multi level, on theta dimension (since um output file size is limited)
_file_types_waroona_ = {'slv':['specific_humidity',  # kg/kg [t,lat,lon]
                               'surface_air_pressure', # Pa 
                               'air_pressure_at_sea_level', # Pa [t,lat,lon]
                               #'surface_altitude', # [lat,lon] # topog # only occurs in first output
                               'surface_temperature', # [t,y,x]
                               'time', # = 6 ;
                               #'level_height', # = 140 ;
                               'latitude', # = 576 ;
                               'longitude', # = 576
                               # x_wind and y_wind are here also...
                               ],
                        'mdl_ro1':['air_pressure', # [t,z,lat,lon]
                               'x_wind', # [t,z,y,x]
                               'y_wind', # [t,z,y,x]
                               #'height_above_reference_ellipsoid', #[z,y,x] # new level heights? # only occurs in first outfile
                               'time', # = 6 ;
                               'level_height', # = 140 ;
                               'latitude', # = 576 ;
                               'longitude', # = 576
                               ],
                        'mdl_th1':['air_pressure', # Pa [t,z,y,x]
                                   'air_temperature', # K [t,z,y,x]
                                   #'height_above_reference_ellipsoid', # m [t,z,y,x] # only occurs in first outfile
                                   'specific_humidity', # kg/kg [tzyx]
                                   'upward_air_velocity', # m/s [tzyx]
                                   'time', # = 6 ;
                                   'level_height', # = 140 ;
                                   'latitude', # = 576 ;
                                   'longitude', # = 576
                                   ],
                        'mdl_th2':['mass_fraction_of_cloud_ice_in_air', # kg/kg [tzyx]
                                   'mass_fraction_of_cloud_liquid_water_in_air',
                                   'time', # = 6 ;
                                   'level_height', # = 140 ;
                                   'latitude', # = 576 ;
                                   'longitude', # = 576
                                   ],
                       }
_files_waroona_ = 'data/waroona/umnsaa_%s_%s.nc' # umnsaa_YYYYMMDDhh_slx.nc


def read_nc_iris(fpath, constraints=None, keepvars=None):
    '''
    Read netcdf file using iris, returning cubeslist
    actual data is not read until you call cube.data
    constraints can be applied, or variable names, or both
    '''
    
    print("INFO: Reading(iris) ",fpath)
    # First read all cubes in file using constraints
    if constraints is not None:
        cubes=iris.load(fpath, constraints)
    else:
        cubes = iris.load(fpath)
    
    # If just some variables are wanted, pull them out
    if keepvars is not None:
        cubes = cubes.extract(keepvars)
    
    return cubes
    
def read_fire(dtimes, constraints=None, extent=None, firefront=True, sensibleheat=False):
    '''
    Read fire output cubes matching dtimes time dim
    '''
    d0 = dtimes[0]
    if d0 < datetime(2016,1,6,15):
        ffpath = 'data/waroona_fire/firefront.CSIRO_24h.20160105T1500Z.nc'
        fluxpath = 'data/waroona_fire/sensible_heat.CSIRO_24h.20160105T1500Z.nc'
    elif d0 < datetime(2016,1,7,15):
        ffpath = 'data/waroona_fire/firefront.CSIRO_24h.20160106T1500Z.nc'
        fluxpath = 'data/waroona_fire/sensible_heat.CSIRO_24h.20160106T1500Z.nc'
    else:
        ffpath = 'data/waroona_fire/firefront.CSIRO_24h.20160107T1500Z.nc'
        fluxpath = 'data/waroona_fire/sensible_heat.CSIRO_24h.20160107T1500Z.nc'
    
    if extent is not None:
        West,East,South,North = extent
        constr_lons = iris.Constraint(longitude = lambda cell: West <= cell <= East )
        constr_lats = iris.Constraint(latitude = lambda cell: South <= cell <= North )
        if constraints is not None:
            constraints = constraints & constr_lats & constr_lons
        else:
            constraints = constr_lats & constr_lons
    
    # Build up cubelist based on which files you want to read
    cubelist = iris.cube.CubeList()
    if firefront:
        ff, = read_nc_iris(ffpath, constraints=constraints)
        ff = subset_time_iris(ff,dtimes,seconds=True)
        cubelist.append(ff)
    if sensibleheat:
        shf, = read_nc_iris(fluxpath, constraints=constraints)
        shf = subset_time_iris(shf,dtimes,seconds=True)
        cubelist.append(shf)
    
    return cubelist
    
def read_z(fpath='data/waroona/umnsaa_2016010515_mdl_th1.nc'):
    '''
    '''
    # read the cube
    z, = read_nc_iris(fpath,'height_above_reference_ellipsoid')
    #height_above_reference_ellipsoid / (m) (model_level_number: 140; latitude: 576; longitude: 576)
    return z
    
def read_waroona(dtime, constraints=None, extent=None, 
                 add_winds=False, add_theta=False, add_dewpoint=True):
    '''
        Read the converted waroona model output files
        returns list of 4 iris cube lists:
        ========
        0: specific_humidity / (1)             (time: 6; latitude: 576; longitude: 576)
        1: surface_air_pressure / (Pa)         (time: 6; latitude: 576; longitude: 576)
        2: air_pressure_at_sea_level / (Pa)    (time: 6; latitude: 576; longitude: 576)
        3: surface_temperature / (K)           (time: 6; latitude: 576; longitude: 576)
        4: topog / (m)                      # added by jwg
        ========
        0: air_pressure / (Pa)                 (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        1: x_wind / (m s-1)                    (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        2: y_wind / (m s-1)                    (time: 6; model_level_number: 140; latitude: 577; longitude: 576)
        3: z_ro / (m)                       # added by jwg
        [ u, v, s / (m s-1) ]               # added here on destaggered horizontal grid if add_winds argument is True
        ========
        0: air_pressure / (Pa)                 (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        1: air_temperature / (K)               (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        2: specific_humidity / (kg kg-1)       (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        3: upward_air_velocity / (m s-1)       (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        4: z_th / (m)                       # added by jwg
        [ potential_temperature / (K) ]     # added if add_theta argument is True
        [ dewpoint_temperature / (K) ]      # added if add_dewpoint is True
        ========
        0: mass_fraction_of_cloud_ice_in_air / (kg kg-1) (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        1: mass_fraction_of_cloud_liquid_water_in_air / (kg kg-1) (time: 6; model_level_number: 140; latitude: 576; longitude: 576)
        2: qc / (g kg-1)                    # added by jwg
    '''
    dstamp=dtime.strftime('%Y%m%d%H')
    
    # First read topography
    topog, = read_nc_iris('data/waroona/umnsaa_2016010515_slv.nc',
                          constraints = 'surface_altitude')
    # also read 3d model height
    zro, = iris.load('data/waroona/umnsaa_2016010515_mdl_ro1.nc', ['height_above_reference_ellipsoid'])
    zth, = iris.load('data/waroona/umnsaa_2016010515_mdl_th1.nc', ['height_above_reference_ellipsoid'])
    
    # If we just want a particular extent, subset to that extent using constraints
    if extent is not None:
        West,East,South,North = extent
        constr_lons = iris.Constraint(longitude = lambda cell: West <= cell <= East )
        constr_lats = iris.Constraint(latitude = lambda cell: South <= cell <= North )
        topog = topog.extract(constr_lons & constr_lats) # subset topog
        if constraints is not None:
            constraints = constraints & constr_lats & constr_lons
        else:
            constraints = constr_lats & constr_lons
    
    if constraints is not None:
        zro = zro.extract(constraints)
        zth = zth.extract(constraints)
    
    cubelists = []
    for filetype,varnames in _file_types_waroona_.items():
        path = _files_waroona_%(dstamp,filetype)
        cubelists.append(read_nc_iris(path,constraints=constraints,keepvars=varnames))
    
    # Add cloud parameter at g/kg scale
    water, ice = cubelists[3].extract(['mass_fraction_of_cloud_liquid_water_in_air',
                                            'mass_fraction_of_cloud_ice_in_air'])
    qc = (water+ice) * 1000
    qc.units = 'g kg-1'
    qc.var_name = 'qc'
    cubelists[3].append(qc)
    # add topog cube
    iris.std_names.STD_NAMES['topog'] = {'canonical_units': 'm'}
    topog.standard_name = 'topog'
    cubelists[0].append(topog)
    # add zro and zth cubes
    iris.std_names.STD_NAMES['z_ro'] = {'canonical_units': 'm'}
    iris.std_names.STD_NAMES['z_th'] = {'canonical_units': 'm'}
    zth.standard_name = 'z_th'
    zro.standard_name = 'z_ro'
    cubelists[1].append(zro)
    cubelists[2].append(zth)
    
    if add_dewpoint:
        # Take pressure and relative humidity
        
        p,q = cubelists[2].extract(['air_pressure','specific_humidity'])
        p_orig_units = p.units
        q_orig_units = q.units
        p.convert_units('hPa')
        q.convert_units('kg kg-1')
        #print("DEBUG: units", p_orig_units, p.units, q_orig_units,q.units)
        # calculate vapour pressure:
        epsilon = 0.6220 # gas constant ratio for dry air to water vapour
        e = p*q/(epsilon+(1-epsilon)*q)
        iris.std_names.STD_NAMES['vapour_pressure'] = {'canonical_units': 'hPa'}
        e.standard_name='vapour_pressure'
        cubelists[2].append(e)
        # calculate dewpoint from vapour pressure
        Td = 234.5 / ((17.67/np.log(e.data/6.112))-1) # in celcius
        Td = Td + 273.16 # celcius to kelvin
        
        # change Td to a Cube
        iris.std_names.STD_NAMES['dewpoint_temperature'] = {'canonical_units': 'K'}
        cubeTd = iris.cube.Cube(Td, standard_name="dewpoint_temperature", 
                                   var_name="Td", units="K",
                                   dim_coords_and_dims=[(p.coord('time'),0),
                                                        (p.coord('model_level_number'),1),
                                                        (p.coord('latitude'),2),
                                                        (p.coord('longitude'),3)])

        cubelists[2].append(cubeTd)
        
    
    if add_winds:
        # wind speeds need to be interpolated onto non-staggered latlons
        p, u1, v1 = cubelists[1].extract(['air_pressure','x_wind','y_wind'])
        
        ### DESTAGGER u and v using iris interpolate
        ### (this will trigger the delayed read)
        # u1: [t,z, lat, lon1]
        # v1: [t,z, lat1, lon]  # put these both onto [t,z,lat,lon]
        u = u1.interpolate([('longitude',p.coord('longitude').points)],
                           iris.analysis.Linear())
        v = v1.interpolate([('latitude',p.coord('latitude').points)],
                           iris.analysis.Linear())
        # Get wind speed cube using hypotenuse of u,v
        s = iris.analysis.maths.apply_ufunc(np.hypot,u,v)
        s.units = 'm s-1'
        # add standard names for these altered variables:
        iris.std_names.STD_NAMES['u'] = {'canonical_units': 'm s-1'}
        iris.std_names.STD_NAMES['v'] = {'canonical_units': 'm s-1'}
        u.standard_name='u'
        v.standard_name='v'
        s.var_name='s' # s doesn't come from a var with a std name so can just use var_name
        cubelists[1].append(u)
        cubelists[1].append(v)
        cubelists[1].append(s)
    
    if add_theta:
        # Estimate potential temp
        p, Ta = cubelists[2].extract(['air_pressure','air_temperature'])
        theta = potential_temperature(p.data,Ta.data)
        # create cube 
        iris.std_names.STD_NAMES['potential_temperature'] = {'canonical_units': 'K'}
        cubetheta = iris.cube.Cube(theta, standard_name="potential_temperature", 
                                   var_name="theta", units="K",
                                   dim_coords_and_dims=[(p.coord('time'),0),
                                                        (p.coord('model_level_number'),1),
                                                        (p.coord('latitude'),2),
                                                        (p.coord('longitude'),3)])
        cubelists[2].append(cubetheta)
    return cubelists

def read_waroona_pcfile(dtime, constraints=None, extent=None, add_winds=False, add_theta=False, add_dewpoint=False):
    '''
    WARNING: add_dewpoint changes air_pressure units to hPa
    0: mass_fraction_of_cloud_liquid_water_in_air / (kg kg-1) (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    1: mass_fraction_of_cloud_ice_in_air / (kg kg-1) (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    2: upward_air_velocity / (m s-1)       (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    3: air_pressure / (hPa)                (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    4: air_temperature / (K)               (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    5: x_wind / (m s-1)                    (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    6: y_wind / (m s-1)                    (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    7: air_pressure_at_sea_level / (Pa)    (time: 2; latitude: 88; longitude: 106)
    8: specific_humidity / (kg kg-1)       (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    9: qc / (g kg-1)                       (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    10: topog / (m)                         (latitude: 88; longitude: 106)
    11: z_th / (m)                          (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    ## Added if add_winds is True
        12: u / (m s-1)                    (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
        13: v / (m s-1)                    (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
        14: s / (m s-1)                    (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    ## Added if add_dewpoint is True
        15: vapour_pressure / (100 Pa)     (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
        16: dewpoint_temperature / (K)     (time: 2; model_level_number: 140; latitude: 88; longitude: 106)
    '''
    dstamp=dtime.strftime('%Y%m%d%H')
     
    # First read topography
    topog, = read_nc_iris('data/waroona/umnsaa_pa2016010515.nc',
                          constraints = 'surface_altitude')
    
    
    # If we just want a particular extent, subset to that extent using constraints
    if extent is not None:
        West,East,South,North = extent
        constr_lons = iris.Constraint(longitude = lambda cell: West <= cell <= East )
        constr_lats = iris.Constraint(latitude = lambda cell: South <= cell <= North )
        topog = topog.extract(constr_lons & constr_lats) # subset topog
        if constraints is not None:
            constraints = constraints & constr_lats & constr_lons
        else:
            constraints = constr_lats & constr_lons
    
    path = 'data/waroona/umnsaa_pc%s.nc'%dstamp
    varnames = ['mass_fraction_of_cloud_liquid_water_in_air',
                'mass_fraction_of_cloud_ice_in_air',
                'upward_air_velocity', # m/s [t,z,lat,lon]
                'air_pressure', # Pa [t,z,lat,lon]
                'air_temperature', # K [lat,lon]
                'air_temperature_0', # K [t, z, lat, lon]
                'x_wind','y_wind', # m/s [t,z,lat,lon]
                'air_pressure_at_sea_level', # Pa [time, lat, lon]
                'specific_humidity', # kg/kg [lat,lon]
                'specific_humidity_0', # kg/kg [t,z,lat,lon] #???
                ]
    cubes = read_nc_iris(path,constraints=constraints,keepvars=varnames)
    
    # specific_humidity is the name of two variables, we just want the good one
    sh0,sh = cubes.extract('specific_humidity')
    if len(sh.shape)==2:
        cubes.remove(sh)
    else:
        cubes.remove(sh0)
    # air_temperature is the name of two variables, we just want the good one
    Ta0,Ta = cubes.extract('air_temperature')
    if len(Ta.shape)==2:
        cubes.remove(Ta)
    else:
        cubes.remove(Ta0)
    
    # Add cloud parameter at g/kg scale
    water, ice = cubes.extract(['mass_fraction_of_cloud_liquid_water_in_air',
                                'mass_fraction_of_cloud_ice_in_air'])
    qc = (water+ice) * 1000
    qc.units = 'g kg-1'
    qc.var_name = 'qc'
    cubes.append(qc)
    # add topog cube
    iris.std_names.STD_NAMES['topog'] = {'canonical_units': 'm'}
    topog.standard_name = 'topog'
    cubes.append(topog)
    # add zth cube
    p,pmsl = cubes.extract(['air_pressure','air_pressure_at_sea_level'])
    nt,nz,ny,nx = p.shape
    reppmsl = np.repeat(pmsl.data[:,np.newaxis,:,:],nz, axis=1) # repeat surface pressure along z axis
    zth = -(287*300/9.8)*np.log(p.data/reppmsl)
    iris.std_names.STD_NAMES['z_th'] = {'canonical_units': 'm'}
    zthcube=iris.cube.Cube(zth, standard_name='z_th', 
                           var_name="zth", units="m",
                           dim_coords_and_dims=[(p.coord('time'),0),
                                                (p.coord('model_level_number'),1),
                                                (p.coord('latitude'),2),
                                                (p.coord('longitude'),3)])
    cubes.append(zthcube)
    
    if add_winds:
        # wind speeds need to be interpolated onto non-staggered latlons
        p, u1, v1 = cubes.extract(['air_pressure','x_wind','y_wind'])
        
        ### DESTAGGER u and v using iris interpolate
        ### (this will trigger the delayed read)
        # u1: [t,z, lat, lon1]
        # v1: [t,z, lat1, lon]  # put these both onto [t,z,lat,lon]
        u = u1.interpolate([('longitude',p.coord('longitude').points)],
                           iris.analysis.Linear())
        v = v1.interpolate([('latitude',p.coord('latitude').points)],
                           iris.analysis.Linear())
        # Get wind speed cube using hypotenuse of u,v
        s = iris.analysis.maths.apply_ufunc(np.hypot,u,v)
        s.units = 'm s-1'
        # add standard names for these altered variables:
        iris.std_names.STD_NAMES['u'] = {'canonical_units': 'm s-1'}
        iris.std_names.STD_NAMES['v'] = {'canonical_units': 'm s-1'}
        u.standard_name='u'
        v.standard_name='v'
        s.var_name='s' # s doesn't come from a var with a std name so can just use var_name
        cubes.append(u)
        cubes.append(v)
        cubes.append(s)
    
    if add_dewpoint:
        # Take pressure and relative humidity
        
        p,q = cubes.extract(['air_pressure','specific_humidity'])
        p_orig_units = p.units
        q_orig_units = q.units
        p.convert_units('hPa')
        q.convert_units('kg kg-1')
        #print("DEBUG: units", p_orig_units, p.units, q_orig_units,q.units)
        # calculate vapour pressure:
        epsilon = 0.6220 # gas constant ratio for dry air to water vapour
        e = p*q/(epsilon+(1-epsilon)*q)
        iris.std_names.STD_NAMES['vapour_pressure'] = {'canonical_units': 'hPa'}
        e.standard_name='vapour_pressure'
        cubes.append(e)
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

        cubes.append(cubeTd)
    
    if add_theta:
        # Estimate potential temp
        p, Ta = cubes.extract(['air_pressure','air_temperature'])
        theta = potential_temperature(p.data,Ta.data)
        # create cube 
        iris.std_names.STD_NAMES['potential_temperature'] = {'canonical_units': 'K'}
        cubetheta = iris.cube.Cube(theta, standard_name="potential_temperature", 
                                   var_name="theta", units="K",
                                   dim_coords_and_dims=[(p.coord('time'),0),
                                                        (p.coord('model_level_number'),1),
                                                        (p.coord('latitude'),2),
                                                        (p.coord('longitude'),3)])
        cubes.append(cubetheta)
    return cubes
    
def read_waroona_oldold(constraints=None, extent=None,
                        add_winds=False):
    '''
    '''
    ddir = model_outputs['waroona_oldold']['path']
    # TODO update to NCI path
    xwind_path = ddir+'oldold_xwind_s5_subset.nc'
    ywind_path = ddir+'oldold_ywind_s5_subset.nc'
    zwind_path = ddir+'oldold_zwind_s5_subset.nc'
    
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
        West,East,South,North = extent
        constr_lons = iris.Constraint(longitude = lambda cell: West <= cell <= East )
        constr_lats = iris.Constraint(latitude = lambda cell: South <= cell <= North )
    
        if constraints is None:    
            constraints = constr_lats & constr_lons
        else:
            constraints = constraints & constr_lats & constr_lons
    
        
    xwind1, = read_nc_iris(xwind_path)#,constraints=constraints)
    ywind1, = read_nc_iris(ywind_path)#,constraints=constraints)
    zwind, = read_nc_iris(zwind_path)#,constraints=constraints)
    
    topog1, = read_nc_iris(ddir+model_outputs['waroona_oldold']['topog'])
    topog = topog1[0,0] # length one time dim and levels dim removed
    
    
    ## Read heights of theta and rho levels
    zth1, = read_nc_iris(ddir+'stage5_ml_htheta.nc')
    zth   = zth1[0] # remove length 1 time dim
    zrho1, = read_nc_iris(ddir+'stage5_ml_hrho.nc')
    zrho   = zrho1[0] # remove length 1 time dim
    iris.std_names.STD_NAMES['z_th'] = {'canonical_units': 'm'}
    iris.std_names.STD_NAMES['z_rho'] = {'canonical_units': 'm'}
    zth.standard_name='z_th'
    zrho.standard_name='z_rho'
    
    # add latlon dims to 2d data
    topog.add_dim_coord(latdim,0)
    topog.add_dim_coord(londim,1)
    
    # add lat/lon dims to data
    xwind1.add_dim_coord(latdim,2) # still staggered
    xwind1.add_dim_coord(londim1,3)
    
    ywind1.add_dim_coord(latdim1,2) # still staggered
    ywind1.add_dim_coord(londim,3)
    
    zwind.add_dim_coord(latdim,2)
    zwind.add_dim_coord(londim,3)
    
    # add lat and lon dims for 3d data
    for cube in [zth, zrho]:
        cube.add_dim_coord(latdim,1)
        cube.add_dim_coord(londim,2)
    
    # add cube to list (subset if possible)
    for cube in [topog, xwind1, ywind1, zwind, zth, zrho]:
        if constraints is not None:
            cube=cube.extract(constraints)
        cubelist.append(cube)
    topog, xwind1, ywind1, zwind, zth, zrho = cubelist # for further calculation
    # update dims if we are subsetting
    if constraints is not None:
        londim = topog.coord('longitude')
        latdim = topog.coord('latitude')
        
    if add_winds:
        
        ### DESTAGGER u and v using iris interpolate
        ### (this will trigger the delayed read)
        # u1: [t,z, lat, lon1]
        # v1: [t,z, lat1, lon]  # put these both onto [t,z,lat,lon]
        u = xwind1.interpolate([('longitude',londim.points)],
                           iris.analysis.Linear())
        v = ywind1.interpolate([('latitude',latdim.points)],
                           iris.analysis.Linear())
        
        # Get wind speed cube using hypotenuse of u,v
        s = iris.analysis.maths.apply_ufunc(np.hypot,u,v)
        s.units = 'm s-1'
        # add standard names for these altered variables:
        iris.std_names.STD_NAMES['u'] = {'canonical_units': 'm s-1'}
        iris.std_names.STD_NAMES['v'] = {'canonical_units': 'm s-1'}
        u.standard_name='u'
        v.standard_name='v'
        s.var_name='s' # s doesn't come from a var with a std name so can just use var_name
        
        cubelist.append(u)
        cubelist.append(v)
        cubelist.append(s)
    
    return cubelist

def read_topog(extentname='waroona'):
    '''
    Read topography cube
    '''
    fpath = 'data/waroona/umnsaa_2016010515_slv.nc'
    if extentname=='sirivan':
        fpath='data/sirivan/umnsaa_pa2017021121.nc'
    # First read topography
    topog, = read_nc_iris(fpath,
                          constraints = 'surface_altitude')
    
    return topog

def subset_time_iris(cube,dtimes,seconds=True,seccheck=121):
    '''
    take a cube with the time dimension and subset it to just have dtimes
    can handle iris seconds or hours time dim formats
    assert times are available within seccheck seconds of dtimes
    '''
    tdim  = cube.coord('time')
    secmult = 1
    if seconds:
        unitformat = 'seconds since %Y-%m-%d %H:%M:00'
    else:
        unitformat = 'hours since %Y-%m-%d %H:%M:00'
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
    assert np.all([(dt[tinds][i] - dtimes[i]).seconds < seccheck for i in range(len(dtimes))]), "fire times are > 2 minutes from requested dtimes"
    
    # subset cube to desired times
    return cube[tinds]

#======================================================
#======================== SIR IVAN STUFF ==============
#======================================================

def read_pcfile(fpath, keepvars=None):
    '''
    Read a umnsaa_pc file, add some extra things such as potential temperature
    '''
    variables=read_nc(fpath,keepvars)
    
    if 'air_pressure' in variables.keys():
        p   = variables['air_pressure'][:,:,:,:] # [time,lev,lat,lon] in pascals
        nt,nz,ny,nx = p.shape

    # height data (z) based on P = P0 exp{ -z/H } with assumed scale height H
    if np.all( [varname in variables.keys() for varname in ['air_pressure','air_pressure_at_sea_level']] ):
        pmsl = variables['air_pressure_at_sea_level'][0,:,:] # [lev,lat,lon]
        zth = np.zeros(np.shape(p))
        for tstep in range(2):
            zth[tstep] = -(287*300/9.8)*np.log(p[tstep]/pmsl[np.newaxis,:,:])
    
        variables['zth']= zth

    # Potential temperature based on https://en.wikipedia.org/wiki/Potential_temperature
    # with gas constant R = 287.05 and specific heat capacity c_p = 1004.64
    if np.all( [varname in variables.keys() for varname in ['air_pressure','air_temperature']] ):
        theta = np.zeros(np.shape(p))
        Ta  = variables['air_temperature'][:,:] # [lat,lon] at surface
        repTa = np.repeat(Ta[np.newaxis,:,:], nz, axis=0) # repeat Ta along z axis
        for tstep in range(2):
            theta[tstep] = repTa*(1e5/p[tstep])**(287.05/1004.64)
        variables['theta'] = theta

    ## Destagger winds
    #u = np.tile(np.nan,(nz,ny,nx))
    if  np.all( [varname in variables.keys() for varname in ['x_wind','y_wind']] ):
        u1  = variables['x_wind'][:,:,:,:] #[time,levs,lat,lon1] wind speeds are on their directional grid edges
        v1  = variables['y_wind'][:,:,:,:] #[time,levs,lat1,lons]
        u = np.tile(np.nan,u1.shape) # tile repeats the nan accross nz,ny,nx dimensions
        u[:,:,:,1:] = 0.5*(u1[:,:,:,1:] + u1[:,:,:,:-1]) # interpolation of edges
        v = 0.5*(v1[:,:,1:,:] + v1[:,:,:-1,:]) # interpolation of edges
        s = np.hypot(u,v) # Speed is hypotenuse of u and v
        ## HYPOT on this 
        # S[0,:,:,0] ARE ALL -5000
        # S[1,:,:,0] ARE ALL NaN
        s[:,:,:,0] = np.NaN # set that edge to NaN
    
        variables['x_wind_destaggered'] = u
        variables['y_wind_destaggered'] = v
        variables['wind_speed'] = s

    
    return variables

def read_sirivan(fpaths, toplev=80, keepvars=None):
    '''
    Read a subset of the Sir Ivan model outputs
    This will loop over the desired files, so RAM doesn't need to be more than 8GB hopefully
    
    INPUTS:
        fpaths is a list of files to read
        toplev = highest model level to read (save some ram)
        keepvars
            set this if you just one one or two outputs over time
    '''
    if isinstance(fpaths, str):
        fpaths = [fpaths]
    # dictionary of data to keep
    data = {}
    
    ## first grab topography from pa vile
    topog,latt,lont = read_topog(_topog_sirivan_)
    data['topog']=topog
    
    ## read first file
    data0 = read_pcfile(fpaths[0],keepvars=keepvars)
    
    # Keep dimensions:
    dims = ['latitude','longitude',
            'time_1', # how many time steps in file 
            'time_0', # time step of file
            'model_level_number', # 140 levels
            'pseudo_level', # 5 levels for soil?
            'longitude_0','latitude_0', # staggered lat,lon
            'time_bnds', # file time bounds
            ]
    # copy dims
    for dim in dims:
        if dim in data0.keys():
            data[dim] = data0[dim]
    
    #assert np.all(data['latitude'] == latt), "topog latitude doesn't match pc file latitude"
    #assert np.all(data['longitude'] == lont), "topog longitude doesn't match pc file longitude"
    if (not np.all(data['latitude'] == latt)) or (not np.all(data['longitude'] == lont)):
        data['latt'] = latt
        data['lont'] = lont
    
    ## data without time dim
    flatdata = ['air_temperature',      # [lat, lon] K at surface
                'specific_humidity',    # [lat,lon] kg/kg
                'surface_air_pressure', # [lat,lon] Pa
                'surface_temperature',  # [lat,lon] K
                #'level_height',         # [z] metres from ground level]
                ]
    if keepvars is not None:
        flatdata = list(set(flatdata) & set(keepvars))
    
    ## 4 dim data [t1, z, lat, lon]
    fulldata = ['air_pressure',                                 # Pascals
                'air_temperature_0',                            # K
                'mass_fraction_of_cloud_ice_in_air',            # kg/kg
                'mass_fraction_of_cloud_liquid_water_in_air',   # kg/kg
                'theta',                                        # K (potential temp)
                'wind_speed',                                   # m/s (wind speed)
                'specific_humidity_0',                          # kg/kg
                'upward_air_velocity',                          # m/s
                'x_wind_destaggered',                           # m/s
                'y_wind_destaggered',                           # m/s
                'zth',                                          # m (height approximated)
                ]
    if keepvars is not None:
        fulldata = list(set(fulldata) & set(keepvars))
    
    ## Start tracking the time dimensions
    hour = []
    time = []
    
    # Gregorian hours since 19700101 00:00:00
    if 'time_0' in data0.keys():
        hour.append(data0['time_0']) 
    else:
        hour.append(np.nanmean(data0['time_1']))
    time.append(data0['time_1'][0])
    time.append(data0['time_1'][1])
    
    ## Copy all the data we want
    for varname in flatdata:
        # add time dim to the flat data
        data[varname] = data0[varname][np.newaxis]
    ## And read the stuff with time dim normally
    for varname in fulldata:
        data[varname] = data0[varname]
        
    ## also copy height
    if 'level_height' in data0.keys():
        data['level_height'] = data0['level_height'][:toplev]
    
    del data0
    
    # read each file and append along time dimension
    for i,fpath in enumerate(fpaths):
        if i==0:
            continue
        datai = read_pcfile(fpath,keepvars=keepvars)
        
        for varname in flatdata:
            # add time dim to the flat data
            datai[varname] = datai[varname][np.newaxis]
        
        # Append each file along the time dimension, which has been added to the vars missing a time dim
        for varname in flatdata+fulldata:
            print("DEBUG:",varname,data[varname].shape, datai[varname].shape)
            data[varname] = np.append(data[varname], datai[varname], axis=0)
            print("DEBUG:",varname,data[varname].shape)
        
        # Append time step
        hour.append(datai['time_0']) 
        time.append(datai['time_1'][0])
        time.append(datai['time_1'][1])
        
        del datai


    
    data['hour'] = np.array(hour)
    data['time'] = np.array(time)
    # also for convenience add cloud metric
    if np.all( [varname in data.keys() for varname in ['mass_fraction_of_cloud_ice_in_air', 'mass_fraction_of_cloud_liquid_water_in_air']] ):
        data['qc'] = data['mass_fraction_of_cloud_ice_in_air'] + data['mass_fraction_of_cloud_liquid_water_in_air']


    # Finally chop the top
    for varname in fulldata:
		# Just read up to toplev 
        data[varname] = data[varname][:,:toplev,:,:]
    

    return data