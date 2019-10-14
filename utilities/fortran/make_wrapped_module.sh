#!/bin/bash

## First load modules which enable python to wrap the fortran code
#module use /g/data3/hh5/public/modules
#module load conda/analysis3

## need fortran compilers? they may be part of python environment
module load python3

# originally I ran python -m numpy.f2py -c HFDiag_calc.f90 -m HFcalc
# now just running 
f2py3 -c HFDiag_calc.f90 -m HFcalc

echo "maybe you need to rename the created file to HFcalc.so"
echo "This will allow you to import HFcalc"
