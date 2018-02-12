#!/bin/sh

#PBS -N TVBii
#PBS -m abe
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -l vmem=30GB

###########################################################################
# Initiate TVBii script to run on HPC infrastructure.
# 
# Written by Hannelore Aerts 
# (UGent, Faculty of Psychology, Department of Data Analysis)
# Date last modification: 28/11/2017
#
# Instructions: type in HPC terminal "qsub -t 1-n HPC_start_version.sh"
###########################################################################

# Set directories
export project="DirName"
export HOMEDIR=~/scratch_vo_user/$project
cd $HOMEDIR

# Run baby, run!
./TVBii_run.sh ${PBS_ARRAYID}

echo "Jobs finished"


