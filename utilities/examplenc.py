#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 15:12:31 2019

@author: jesse
"""

from netCDF4 import Dataset
import numpy as np


####
#### Example writing a netcdf file (https://unidata.github.io/netcdf4-python/netCDF4/index.html)
####
# open for write
rootgrp = Dataset("data/test.nc", "w", format="NETCDF4")

# create metadata?

# create some groups
fcstgrp = rootgrp.createGroup("forecasts")
analgrp = rootgrp.createGroup("analyses")

fcstgrp1 = rootgrp.createGroup("/forecasts/model1")
fcstgrp2 = rootgrp.createGroup("/forecasts/model2")

# show summary
print(rootgrp.groups)

'''
Create dimensions

A Python string is used to set the name of the dimension, and an integer value is used to set the size. To create an unlimited dimension (a dimension that can be appended to), the size value is set to None or 0. In this example, there both the time and level dimensions are unlimited. Having more than one unlimited dimension is a new netCDF 4 feature, in netCDF 3 files there may be only one, and it must be the first (leftmost) dimension of the variable
'''
level = rootgrp.createDimension("level", None)
time = rootgrp.createDimension("time", None)
lat = rootgrp.createDimension("lat", 73)
lon = rootgrp.createDimension("lon", 144)
# stored as dictionary
print(rootgrp.dimensions)
# length shown by len, items have method isunlimited 
print(len(lon),lon.isunlimited())

'''
Create variables

The createVariable method has two mandatory arguments, the variable name (a Python string), and the variable datatype. The variable's dimensions are given by a tuple containing the dimension names (defined previously with createDimension). To create a scalar variable, simply leave out the dimensions keyword.
'''
# Dimensions are also variables
times = rootgrp.createVariable("time","f8",("time",))
levels = rootgrp.createVariable("level","i4",("level",))
latitudes = rootgrp.createVariable("lat","f4",("lat",))
longitudes = rootgrp.createVariable("lon","f4",("lon",))
# two dimensions unlimited
temp = rootgrp.createVariable("temp","f4",("time","level","lat","lon",))
# show summary
print (temp)

# Create variable inside heirarchy
ftemp = rootgrp.createVariable("/forecasts/model1/temp","f4",("time","level","lat","lon",))

# can print data for group or variable path:
print(rootgrp["/forecasts/model1"]) # a Group instance

'''
Create attributes

There are two types of attributes in a netCDF file, global and variable. Global attributes provide information about a group, or the entire dataset, as a whole. Variable attributes provide information about one of the variables in a group. Global attributes are set by assigning values to Dataset or Group instance variables. Variable attributes are set by assigning values to Variable instances variables. Attributes can be strings, numbers or sequences. Returning to our example,
'''
import time
rootgrp.description = "bogus example script"
rootgrp.history = "Created " + time.ctime(time.time())
rootgrp.source = "netCDF4 python module tutorial"
latitudes.units = "degrees north"
longitudes.units = "degrees east"
levels.units = "hPa"
temp.units = "K"
times.units = "hours since 2001-01-01 00:00:00.0"
times.calendar = "gregorian"
# print all attributes
for name in rootgrp.ncattrs():
    print("Global attr", name, "=", getattr(rootgrp,name))
    
# Attributes can be deleted from a netCDF Dataset, Group or Variable using the python del statement (i.e. del grp.foo removes the attribute foo the the group grp).  

# Populate the dimensions
lats =  np.arange(-90,91,2.5)
lons =  np.arange(-180,180,2.5)
latitudes[:] = lats
longitudes[:] = lons
print("latitudes =\n",latitudes[:])

# Unlimited dimensions grow with data
nlats = len(rootgrp.dimensions["lat"])
nlons = len(rootgrp.dimensions["lon"])
print("temp shape before adding data = ",temp.shape)
from numpy.random import uniform
temp[0:5,0:10,:,:] = uniform(size=(5,10,nlats,nlons))
print("temp shape after adding data = ",temp.shape)

# levels have grown, but no values yet assigned.
print("levels shape after adding pressure data = ",levels.shape)
# Can do some fancy slicing (not the same as numpy arrays)
tempdat = temp[::2, [1,3,6], lats>0, lons>0]
print("shape of fancy temp slice = ",tempdat.shape)


'''
Time coordinates

Time coordinate values pose a special challenge to netCDF users. Most metadata standards (such as CF) specify that time should be measure relative to a fixed date using a certain calendar, with units specified like hours since YY-MM-DD hh:mm:ss. These units can be awkward to deal with, without a utility to convert the values to and from calendar dates. The function called num2date and date2num are provided with this package to do just that (starting with version 1.4.0, the cftime package must be installed separately). Here's an example of how they can be used:
'''
# fill in times.
from datetime import datetime, timedelta
from netCDF4 import num2date, date2num
dates = [datetime(2001,3,1)+n*timedelta(hours=12) for n in range(temp.shape[0])]
times[:] = date2num(dates,units=times.units,calendar=times.calendar)
print("time values (in units %s): " % times.units+"\n",times[:])

dates = num2date(times[:],units=times.units,calendar=times.calendar)
print("dates corresponding to time values:\n",dates)


# Make sure to flush out to nc file
rootgrp.close()


'''
compound data types.

Compound data types map directly to numpy structured (a.k.a 'record') arrays. Structured arrays are akin to C structs, or derived types in Fortran. They allow for the construction of table-like structures composed of combinations of other data types, including other compound types. Compound types might be useful for representing multiple parameter values at each point on a grid, or at each time and space location for scattered (point) data. You can then access all the information for a point by reading one variable, instead of reading different parameters from different variables. Compound data types are created from the corresponding numpy data type using the createCompoundType method of a Dataset or Group instance. Since there is no native complex data type in netcdf, compound types are handy for storing numpy complex arrays.
'''

'''
Parallel IO.

If MPI parallel enabled versions of netcdf and hdf5 or pnetcdf are detected, and mpi4py is installed, netcdf4-python will be built with parallel IO capabilities enabled. Parallel IO of NETCDF4 or NETCDF4_CLASSIC formatted files is only available if the MPI parallel HDF5 library is available. Parallel IO of classic netcdf-3 file formats is only available if the PnetCDF library is available. To use parallel IO, your program must be running in an MPI environment using mpi4py.
'''
#from mpi4py import MPI
#rank = MPI.COMM_WORLD.rank  # The process ID (integer 0-3 for 4-process run)