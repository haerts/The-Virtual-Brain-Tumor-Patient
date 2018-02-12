#!/bin/bash

#PBS -N TVBii
#PBS -m abe
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -l vmem=30GB

###########################################################################
# Run TVBii code to run on HPC infrastructure.
# 
# Written by Hannelore Aerts 
# (UGent, Faculty of Psychology, Department of Data Analysis)
# Date last modification: 28/11/2017
###########################################################################

if [ -z "$PBS_ARRAYID" ]
then
	echo 'ERROR: $PBS_ARRAYID is not set, submit job(s) using "qsub -t <array expression>"'
	exit 1
fi

# Directory where inputs are located
dos2unix $HOMEDIR/input.txt --quiet

# Make input variable with subIDs
export subID=`sed -n "${PBS_ARRAYID}p" $HOMEDIR/input.txt`
echo "subID: $subID"

# Define TVBii folder
subFolder="$HOMEDIR/subjects/${subID}"
cd $subFolder


# First get TR 
tr=$( cat $subFolder/tr.txt ) 

# Generate parameter files (depending on TR)
LANG=eng_US
it=0

if [ $tr = 2.1 ]; then 
	trm=2100
	simlen=420000
elif [ $tr = 2.4 ]; then 
	trm=2400
	simlen=480000
fi

for G in $(seq 0.01 0.015 3.00)
do
	it=$(($it+1))
	echo "68 $G 0.15 1.40 1.0 0.01 $simlen 10000 $trm 999999.99 42" > param_set_$it
	#parameters: number of nodes, G, J_NMDA, w_plus, J_i, sigma, 
	#simulation length (ms), number of iterations for FIC optimization, BOLD_TR (ms), 
	#global transmission velocity, random seed number
done


# Then copy the TVBii executable
cp $HOMEDIR/tvbii_linux $subFolder/tvbii_linux


# Run TVBii for all parameter files, for all thresholds and distance metrics

# thrA - distLen
for sim in {1..200}
do
	echo "Processing iteration $sim"
	./tvbii_linux param_set_${sim} "${subID}"_scale68_thrA_distLen
done

mv $subFolder/output $subFolder/output_thrA_distLen
mkdir ${subFolder}/output


# thrA - distLog
for sim in {1..200}
do
	echo "Processing iteration $sim"
	./tvbii_linux param_set_${sim} "${subID}"_scale68_thrA_distLog
done

mv $subFolder/output $subFolder/output_thrA_distLog
mkdir ${subFolder}/output


# thrR - distLen
for sim in {1..200}
do
	echo "Processing iteration $sim"
	./tvbii_linux param_set_${sim} "${subID}"_scale68_thrR_distLen
done

mv $subFolder/output $subFolder/output_thrR_distLen
mkdir ${subFolder}/output


# thrR - distLog
for sim in {1..200}
do
	echo "Processing iteration $sim"
	./tvbii_linux param_set_${sim} "${subID}"_scale68_thrR_distLog
done

mv $subFolder/output $subFolder/output_thrR_distLog





