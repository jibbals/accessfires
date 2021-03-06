#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 11:28:52 2020
    Creation and some utilisation of spatially averaged time series metrics
    
@author: jesse
"""
###########
### IMPORTS
###########
import numpy as np
import xarray as xr
import pandas as pd
import iris, os, shutil

from utilities import constants,fio,utils
from fireplan import show_fire_outlines

###########
### GLOBALS
###########
_sn_="metrics"

# altitudes used in metrics file
_level_height_=("level",
                 np.array([0,10, 100, 200, 300, 400, 500, 600, 700, 800, 900, 
                           1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 
                           3500, 4000, 4500, 5000, 6000, 7000, 8000, 9000,
                           10000, 11000, 12000, 13000, 14000]))
_naltitudes_= len(_level_height_[1]) 

###########
### METHODS
###########
def metric_file_path(mr, extentname):
    return "data/metrics/%s_%s.nc"%(mr,extentname)

def metric_file_variables(ntimes=144):
    """
    Set up variables to be used by metric files
    EXAMPLE FILE:
        INFO: reading/plotting  data/metrics/sirivan_run4_sirivanz.nc
        DEBUG: <xarray.Dataset>
        Dimensions:                    (level: 10, pctl: 5, time: 144)
        Coordinates:
          * time                       (time) datetime64[ns] 2017-02-11T11:10:00 ... ...
          * level                      (level) int64 0 1 2 ... 
          * pctl                       (pctl) int64 0 1 2 3 4
        Data variables:
            air_temperature_mean       (time, level) float64 ...
            air_temperature_5ns        (time, level, pctl) float64 ...
            air_pressure_mean          (time, level) float64 ...
            air_pressure_5ns           (time, level, pctl) float64 ...
            wind_direction_mean        (time, level) float64 ...
            wind_direction_5ns         (time, level, pctl) float64 ...
            windspeed_mean             (time, level) float64 ...
            windspeed_5ns              (time, level, pctl) float64 ...
            relative_humidity_mean     (time, level) float64 ...
            relative_humidity_5ns      (time, level, pctl) float64 ...
            level_height               (level) int64 ...
            firespeed_mean             (time) float64 ...
            firespeed_5ns              (time, pctl) float64 ...
            firespeed_nonzero_mean     (time) float64 ...
            firespeed_nonzero_5ns      (time, pctl) float64 ...
            FFDI_mean                  (time) float64 ...
            FFDI_5ns                   (time, pctl) float64 ...
            sensibleheat_mean          (time) float64 ...
            sensibleheat_5ns           (time, pctl) float64 ...
            sensibleheat_nonzero_mean  (time) float64 ...
            sensibleheat_nonzero_5ns   (time, pctl) float64 ...
            firepower                  (time) float64 ...
        Attributes:
            description:  Time series at several model levels (closest index to level...
            WESN:         [149.4  150.19 -32.2  -31.8 ]
    """
    
    dimarraypairs={}
    for varnames in [
            "air_temperature", # temperature
            "air_pressure", # pressure
            "wind_direction", # wind direction
            "windspeed", # wind speed
            "relative_humidity", # rel humidity
            "specific_humidity",
            ]:
        # variable holding mean value
        dimarraypairs[varnames+"_mean"]=(("time","level"),np.zeros([ntimes,_naltitudes_])+np.NaN)
        # variable holding min, Q1, Q2, Q3, max (5 number summary)
        dimarraypairs[varnames+"_5ns"]=(("time","level","pctl"),np.zeros([ntimes,_naltitudes_,5])+np.NaN)
    
    # height coordinate
    dimarraypairs["level_height"]=_level_height_
    
    # some have no level dim:
    for varname in [
            "firespeed", # velocity of fire front? m/s
            "firespeed_nonzero", # with any zeros removed
            "FFDI", # FFDI is single leve, based on 10m winds and surface T,RH
            "sensibleheat",
            "sensibleheat_nonzero",
            ]:
        dimarraypairs[varname+"_mean"]=("time",np.zeros([ntimes])+np.NaN)
        dimarraypairs[varname+"_5ns"]=(("time","pctl"),np.zeros([ntimes,5])+np.NaN)
    dimarraypairs["firepower"]=("time",np.zeros([ntimes])+np.NaN)
    
    return dimarraypairs

def make_empty_metrics_file(mr="sirivan_run4",
                            extentname=None, 
                            **to_netcdf_args):
    """
    create metrics file, can be filled using other method
    ARGS:
        mr: model run name from fio.run_info keys
        extentname: defaults to location+"z" (eg. sirivanz for sirivan_run5)
        other arguments to be sent to xarray.DataSet.to_netcdf, defaults:
            mode="w"
            filename="data/metrics/<mr>_<extentname>.nc"
    
    Data array: variables with _mean and _5ns as follows
        "air_temperature", # temperature (K)
        "air_pressure", # pressure (hPa)
        "wind_direction", # wind direction (meteorological degrees)
        "windspeed", # wind speed (m/s)
        "relative_humidity", # rel humidity (frac?)
        "FFDI", # forest fire danger index
        "firespeed", # fire speed (m/s?)
        "firespeed_nonzero", # zeros removed
        "firepower", # sum of firepower
        Coords:
            "time": time
            "level": model level 0, 10m, 100m, 500m, 1000m, 2000m, 4000m, 6000m, 10000m, 15000m
            "pctl": min, Q1, median, Q3, max
    
    Returns:
        path of created file
    """
    if extentname is None:
        extentname=mr.split("_")[0]+"z"
    
    ## Set up coordinates to be used
    sirivan_10minute_dates = pd.date_range(start="2017-02-11 11:10",end="2017-02-12 11:00",periods=144)
    waroona_10minute_dates = pd.date_range(start="2016-01-05 15:10",end="2016-01-06 15:00",periods=144)
    datecoord = waroona_10minute_dates if "waroona" in mr else sirivan_10minute_dates
    coords = {
        "time": datecoord,
        "level": np.arange(_naltitudes_),
        "pctl":np.arange(5),
    }
    ## Set up variables to be used
    arrays=metric_file_variables(144)
    
    global_attrs={
        "description":"Time series at several model levels (closest index to level_height) created by collapsing horizontally on a subset (bound by [West,East,South,North]) of the model output",
        "WESN":constants.extents[extentname],
        }
    ##  Create dataset using variables and coords
    ds = xr.Dataset(
        data_vars=arrays,
        coords=coords,
        attrs=global_attrs,
        )

    ## Attributes
    key_unit = {"air_temperature":"K",
                "air_pressure":"hPa",
                "windspeed":"m s-1",
                "firespeed":"m s-1",
                "firespeed_nonzero":"m s-1",
                "wind_direction":"degrees",
                "sensibleheat":"Watts m-2",
                "sensibleheat_nonzero":"Watts m-2",
                "specific_humidity":"kg kg-1",
                }
    for key,unit in key_unit.items():
        ds[key+"_mean"].attrs["units"]=unit
        ds[key+"_5ns"].attrs["units"]=unit
    ds["firepower"].attrs["units"]="Gigawatts"
    ds["wind_direction_mean"].attrs["long_name"]="degrees clockwise from due north"
    ds["wind_direction_5ns"].attrs["long_name"]="degrees clockwise from due north"
    ds["FFDI_mean"].attrs["description"]="ffdi=2*np.exp(-0.45 + 0.987*np.log(DF) - .0345*RH+.0338*T+.0234*v): where v is 10 metre wind speed, others are at surface"
    ds["FFDI_5ns"].attrs["description"]="ffdi=2*np.exp(-0.45 + 0.987*np.log(DF) - .0345*RH+.0338*T+.0234*v): where v is 10 metre wind speed, others are at surface"
    
    ## Default arguments to save the netcdf file
    #if "group" not in to_netcdf_args:
    #    to_netcdf_args["group"]=extentname
    if "mode" not in to_netcdf_args:
        to_netcdf_args["mode"]="w"
    if "path" not in to_netcdf_args:
        to_netcdf_args["path"]=metric_file_path(mr,extentname)
    
    fio.make_folder(to_netcdf_args["path"])
    print("INFO: writing ",to_netcdf_args["path"])
    ds.to_netcdf(**to_netcdf_args)
    ds.close()
    return to_netcdf_args["path"]

def make_metrics_from_model(mr,hour=0,extentname=None,HSkip=None):
    """
    return dict of arrays [6, 10, [5]] with _mean and _5ns as follows:
       "air_temperature", # temperature
        "air_pressure", # pressure
        "wind_direction", # wind direction
        "windspeed", # wind speed
        "relative_humidity", # rel humidity
        "FFDI", # forest fire danger index (based on 10m winds and surface T,RH)
        "firepower",
        "firespeed",
    Coords:
        "time": time
        "level": model level 0, 10m, 100m, 500m, 1000m, 2000m, 4000m, 6000m, 10000m, 15000m
        "pctl": min, Q1, median, Q3, max
    """
    if extentname is None:
        extentname=mr.split("_")[0]+"z"
    extent=constants.extents[extentname]
    
    dthour=fio.run_info[mr]["filedates"][hour]
    
    vars_to_make=metric_file_variables(6)
    heights=vars_to_make['level_height'][1]
    
    # Read model data for extentname
    cubes=fio.read_model_run(mr, fdtime=dthour, extent=extent,
                             add_topog=False, HSkip=HSkip)
    # we subset then add wind/rel humidity fields
    #print("DEBUG: read model cubes:",cubes)
    ctimes = utils.dates_from_iris(cubes.extract('air_temperature')[0])
    model_heights = utils.height_from_iris(cubes.extract("air_temperature")[0])
    
    return_dict = {}
    
    ## Get wind speed and direction (after subsetting hopefully)
    u1, v1 = cubes.extract(['x_wind','y_wind']) #staggered but shouldn't matter for timeseries
    #t0 = time.perf_counter()
    u0 = utils.interp_cube_to_altitudes(u1,heights,model_heights, closest=True)
    #t1 = time.perf_counter()
    v0 = utils.interp_cube_to_altitudes(v1,heights,model_heights,closest=True)
    #t2 = time.perf_counter()
    # destagger
    u = u0.interpolate([('longitude',v0.coord('longitude').points)],
                       iris.analysis.Linear())
    v = v0.interpolate([('latitude',u0.coord('latitude').points)],
                       iris.analysis.Linear())
    #print("DEBUG: time with interpolation: %.5f seconds"%(t1-t0))
    #print("DEBUG: time without interpolation: %.5f seconds"%(t2-t1))
    # Get wind speed cube using hypotenuse of u,v
    firesh,firespeed,u10,v10=fio.read_fire(mr, extent=extent, dtimes=ctimes, 
                                             sensibleheat=True, firefront=False, 
                                             wind=True, firespeed=True,
                                             HSkip=HSkip)
    
    s10=utils.wind_speed_from_uv_cubes(u10,v10)
    
    s0 = utils.wind_speed_from_uv_cubes(u,v)
    s=s0.data
    s[:,1,:,:] = s10.data
    cubews = iris.cube.Cube(s,
                            var_name='wind_speed',
                            dim_coords_and_dims=[(s0.coord('time'),0),
                                                 (s0.coord('model_level_number'),1),
                                                 (s0.coord('latitude'),2),
                                                 (s0.coord('longitude'),3)]
                            )
    
    
    # Get wind direction using arctan of y/x
    wd = utils.wind_dir_from_uv(u.data,v.data)
    wd[:,1,:,:] = utils.wind_dir_from_uv(u10.data,v10.data)
    cubewd = iris.cube.Cube(wd,
                            var_name='wind_direction',
                            units='degrees',
                            dim_coords_and_dims=[(s0.coord('time'),0),
                                                 (s0.coord('model_level_number'),1),
                                                 (s0.coord('latitude'),2),
                                                 (s0.coord('longitude'),3)])
    
    # calculate rel humidity
    q1,T1,Pa1 = cubes.extract(['specific_humidity','air_temperature','air_pressure'])
    # compute RH from specific and T in kelvin
    q = utils.interp_cube_to_altitudes(q1,heights,model_heights, closest=True)
    T = utils.interp_cube_to_altitudes(T1,heights,model_heights, closest=True)
    Pa = utils.interp_cube_to_altitudes(Pa1,heights,model_heights, closest=True)
    RH = utils.relative_humidity_from_specific(q.data, T.data)
    iris.std_names.STD_NAMES['relative_humidity'] = {'canonical_units': '1'}
    cubeRH = iris.cube.Cube(RH, standard_name="relative_humidity",
                               var_name="RH", units="1",
                               dim_coords_and_dims=[(q.coord('time'),0),
                                                    (q.coord('model_level_number'),1),
                                                    (q.coord('latitude'),2),
                                                    (q.coord('longitude'),3)])
    
    # also wind speed at 10m for FFDI calc
    WS_10m=np.squeeze(s[:,1,:,:]) # 10m wind speed
    # Surface RH as %
    RH_surf=100*np.squeeze(RH[:,0,:,:])
        
    cubes_subset={
        "wind_direction":cubewd,
        "windspeed":cubews,
        "relative_humidity":cubeRH,
        "firespeed":firespeed,
        "sensibleheat":firesh,
        "specific_humidity":q,
        "air_temperature":T,
        "air_pressure":Pa,
        }
    # also make nonzero metrics from some fire fields
    for cube, varname, units in zip([firespeed,firesh],
                                    ["firespeed","sensibleheat"],
                                    ["m s-1","Watts m-2"]):
        arr_nz = np.copy(cube.data)
        arr_nz[arr_nz<1e-5] = np.NaN
        cube_nz = iris.cube.Cube(arr_nz, 
                                 var_name=varname+"_nonzero", 
                                 units=units,
                                 dim_coords_and_dims=[(firespeed.coord('time'),0),
                                                      (firespeed.coord('latitude'),1),
                                                      (firespeed.coord('longitude'),2)])
        cubes_subset[varname+"_nonzero"]=cube_nz
    
    for varname in [
        "air_temperature", # temperature
        "air_pressure", # pressure
        "wind_direction", # wind direction
        "windspeed", # wind speed
        "relative_humidity", # rel humidity
        "specific_humidity",
        "firespeed", # fire speed
        "firespeed_nonzero",
        "sensibleheat",
        "sensibleheat_nonzero"
        ]:
        print("INFO: collating ",varname)
        
        cube = cubes_subset[varname]
        
        # make contiguous to evade warning
        for coordname in ['latitude','longitude']:
            if not cube.coord(coordname).has_bounds():
                cube.coord(coordname).guess_bounds()
        
        ## horizontal aggregation
        cube_mean = cube.collapsed(['latitude','longitude'], 
                                   iris.analysis.MEAN)
    
        cube_5ns = cube.collapsed(['latitude','longitude'], 
                                  iris.analysis.PERCENTILE, 
                                  percent=[0,25,50,75,100])
        
        return_dict[varname+"_mean"]=cube_mean.data
        # make pctl dimension the last one
        return_dict[varname+"_5ns"]=np.moveaxis(cube_5ns.data, 0, -1)
    
    
    DF=10.0
    print("INFO: collating FFDI (assume DF = %.1f)"%DF)
    # Surface T in Celcius
    T=np.squeeze(cubes.extract('air_temperature')[0][:,0,:,:].data) - 273.15
    # WS in km/h
    FFDI=utils.FFDI(DF,RH_surf,T,WS_10m*3.6)
    
    return_dict["FFDI_mean"]=np.mean(FFDI,axis=(1,2)) # mean over latlon
    return_dict["FFDI_5ns"]=np.moveaxis(np.percentile(FFDI,[0,25,50,75,100],
                                                      axis=(1,2)),
                                        0, -1)
    print("INFO: collating firepower")
    firepower_GW = utils.firepower_from_cube(firesh)
    return_dict["firepower"] = np.sum(firepower_GW,axis=(1,2))
    
    return return_dict
        

def add_to_metrics_file(mr, hour=0, extentname=None, HSkip=2):
    """
    take one hour of model run, get 10 minutely aggregates, update metrics file
    """
    if extentname is None:
        extentname=mr.split("_")[0]+'z'
    fpath = metric_file_path(mr,extentname)
    
    ## May need to create the file
    if not os.path.isfile(fpath):
        fpath = make_empty_metrics_file(mr,extentname)
        
    # temporary path for overwriting purposes
    fpath_tmp = fpath.split(".")[0] + "_tmp.nc"
    
    ## Now we read in the model and create our timeseries
    arrays = make_metrics_from_model(mr,hour=hour,extentname=extentname,HSkip=HSkip)
    nsteps=np.shape(arrays["air_temperature_mean"])[0]
    if nsteps != 6:
        print("WARNING: %d timesteps in output hour!"%nsteps)
    
    # indices for inserting new data
    insertinds=slice(hour*6,hour*6+nsteps)
    
    ## Finally open the dataset, update the metrics, save to temp path, overwrite original path
    with xr.open_dataset(fpath) as ds:
        # load data so we can overwrite it
        ds.load()
        
        # move metrics into file!
        for varname in arrays.keys():
            print("INFO: updating ",varname)
            # time, level, [pctl]
            ds[varname][insertinds] = arrays[varname]
        
        print("INFO: saving to temporary file path ",fpath_tmp)
        ds.to_netcdf(path=fpath_tmp,mode="w")# write to temporary path
    print("INFO: overwriting file path ",fpath)
    shutil.copy(fpath_tmp,fpath) # overwrite original with updated
    os.remove(fpath_tmp) # delete temp file
    
def compare_metrics(mrs=["sirivan_run5",],extent=None,):
    """
    Look at some metrics side by side for a list of model runs
    """
    
    if extent is None:
        extent=mrs[0].split('_')[0]+"z"

    # colours for each model run:
    eight_paired_colours=['#a6cee3','#1f78b4',
                          '#b2df8a','#33a02c',
                          '#fb9a99','#e31a1c',
                          '#fdbf6f','#ff7f00']
    colors = {"sirivan_run3":eight_paired_colours[6],
              "sirivan_run4":eight_paired_colours[7],
              "sirivan_run5":eight_paired_colours[0],
              "sirivan_run5_hr":eight_paired_colours[1],
              "sirivan_run6":eight_paired_colours[2],
              "sirivan_run6_hr":eight_paired_colours[3],
              "sirivan_run7":eight_paired_colours[4],
              "sirivan_run7_hr":eight_paired_colours[5],
              "waroona_run3":eight_paired_colours[3],
              }
    
    
    fig,axes = plt.subplots(5,1,figsize=[12,12])
    
    
    ## loop over runs, adding the data to premade figure axes
    for mr in mrs:
        color=colors[mr]
        fpath=metric_file_path(mr,extent)

        with xr.open_dataset(fpath) as ds:
            print("INFO: reading/plotting ",fpath)
            #print("DEBUG:",ds)
            kwargs={"label":mr,"color":color}
            kwargs_maximums={
                "marker":"^",
                "linestyle":"None",
                "color":color,
                "alpha":0.6,
                "label":'_nolegend_',#no label for maximums
                }
            ## temperature
            plt.sca(axes[0])
            plt.plot(ds["air_temperature_mean"][:,0]-273.15, **kwargs)
            
            ## 10m windspeed
            plt.sca(axes[1])
            plt.plot(ds["s_mean"][:,1], **kwargs)
            plt.plot(ds["s_5ns"][:,1,4], **kwargs_maximums)
            
            plt.sca(axes[2])
            plt.plot(ds["firespeed_nonzero_mean"], **kwargs)
            # max firespeed
            plt.plot(ds["firespeed_nonzero_5ns"][:,4], **kwargs_maximums)
            
            plt.sca(axes[3])
            plt.plot(ds["sensibleheat_nonzero_mean"], **kwargs)
            plt.plot(ds["sensibleheat_nonzero_5ns"][:,4], **kwargs_maximums)
            
            plt.sca(axes[4])
            plt.plot(ds["firepower"], **kwargs)
            
    legflag=False
    for ax,title,ylabel in zip(
            axes,
            ["surface air-temperature", "10m wind speed", "firespeed", 
             "heatflux", "firepower"],
            ["Celcius", "ms$^{-1}$", "ms$^{-1}$", 
             "Watts m$^{-2}$", "Gigawatts"]):
        plt.sca(ax)
        plt.title(title)
        plt.ylabel(ylabel)
        if not legflag:
            plt.legend() 
            legflag=True
    plt.xlabel("time")
    fio.save_fig("comparison",_sn_,"metrics_%s.png"%extent,plt)
    

if __name__=="__main__":
    import matplotlib.pyplot as plt
    
    siruns=["sirivan_run4","sirivan_run5","sirivan_run5_hr",
            "sirivan_run6","sirivan_run6_hr",
            "sirivan_run7","sirivan_run7_hr"]
    wruns=["waroona_run3",]

    ## Check plots
    if False:
        compare_metrics(mrs=wruns)
    
    ## metric file creation/population ~ 90GB RAM and 1:20:00 CPU time
    if True:
        mr = "waroona_run3"
        extent="waroonaz"
        fpath = make_empty_metrics_file(mr,extentname=extent)
        for hour in range(1):
            add_to_metrics_file(mr, hour=hour, extentname=extent, HSkip=2)
    
    ## Show area where metrics are averaged
    if False:
        for mr, extentname in zip(['waroona_run3','sirivan_run6_hr'],['waroonaz','sirivanz']):
            fig,ax = show_fire_outlines(mr, extentname)
            fio.save_fig("comparison", _sn_, extentname+'_'+mr+'_extent.png', plt)
