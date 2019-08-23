#!/bin/bash
#PBS -P en0
#PBS -N loop_making
#PBS -q express
#PBS -l walltime=4:00:00
#PBS -l mem=40GB
#PBS -l wd
#PBS -l ncpus=1
#PBS -l software=python


## Check running as script on queue
if [ -z ${PBS_O_LOGNAME} ] || [ -z ${quarterday} ] 
then
    echo "EG usage: qsub -v quarterday=2 ${0}"
    echo "    to run make_loops.py on 2nd 6 hour period"
    exit 0
fi


# python environment including jupyter and plotting stuff
# module load python3/3.7.2
# public module has conda environment with jupyter, cartopy, and lots of stuff
module use /g/data3/hh5/public/modules
module load conda/analysis3

# run script with winds flag
./make_loops.py --quarterday ${quarterday} --winds
