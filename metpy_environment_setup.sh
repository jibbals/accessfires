#!/bin/bash

# This stuff is what I used to create an environment on NCI similar to what I use on my laptop:
# Laptop environmetn was created using 
#   conda create -n BOMenv python=3.7 spyder matplotlib cartopy numpy scipy
# Running this script with source environment_setup.sh should give you the right modules

# don't need spyder (python GUI) on NCI, just use vim or whatever editor you like

# Python3, includes matplotlib, numpy, and scipy?
#module load python3/3.7.2
# access module has conda environment with jupyter, cartopy, and lots of stuff
module use /g/data3/hh5/public/modules
# Specific version required for using metpy
module load conda/analysis3-19.07

# module for some threedee stuff (warnings are thrown otherwise)
#module load gtk

