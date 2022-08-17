#!/bin/bash

# To run this script, do 
# qsub -t 1-n submit_permutationTest.sh CONFIG IDS PERMUTATION_NUMBER CORES
#
#$ -N circPerm
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=64G
#$ -l h_rt=24:00:00
#$ -pe sharedmem 12

CONFIG=$1
IDS=$2
PERMUTATION_NUMBER=$3
CORES=$4

source $CONFIG

SAMPLE=$(head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 1)
INPUTA=$(head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 2)
INPUTB=$(head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 3)


export PATH=/exports/eddie/scratch/s1749179/anaconda/envs/stats/bin:$PATH

PERMUTATION=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/src/scripts/permutation.R

cd $OVERLAP_STATS

Rscript $PERMUTATION \
	--inputA ${MACS2_PEAKS}/${INPUTA} \
	--inputB ${MACS2_PEAKS}/${INPUTB} \
	--sample $SAMPLE \
	--ntimes $PERMUTATION_NUMBER \
	--cores $CORES
