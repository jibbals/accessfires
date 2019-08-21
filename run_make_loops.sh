#!/bin/bash
#PBS -P en0
#PBS -N loop_making
#PBS -q express
#PBS -l walltime=2:00:00
#PBS -l mem=150GB
#PBS -l wd
#PBS -l ncpus=4
#PBS -l software=python


## Check running as script on queue
if [ -z ${PBS_O_LOGNAME} ] ; then
    echo "EG usage: qsub ${0}"
    echo "    to run make_loops.py"
    exit 0
fi


# python environment including jupyter and plotting stuff
module load python3/3.7.2
# public module has conda environment with jupyter, cartopy, and lots of stuff
module use /g/data3/hh5/public/modules
module load conda/analysis3

# run script and send output to out.make_loops, also send error stream to output stream
python make_loops.py
