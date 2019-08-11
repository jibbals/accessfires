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

#import numpy as np
from netCDF4 import Dataset


def read_nc(fpath, keepvars=None):
  '''
    Generic read function for netcdf files
    Reads all dimensions, and all [or some] variables into a dictionary to be returned
  '''
  ncfile =  Dataset(fpath,'r')
  variables = {}
  if keepvars is None:
    for vname in ncfile.variables:
      variables[vname] = ncfile.variables[vname]
  else:
    # keep dimensions and whatever is in keepvars
    for vname in set(keepvars) | set(ncfile.dimensions.keys()):
      variables[vname] = ncfile.variables[vname][:]
  return variables