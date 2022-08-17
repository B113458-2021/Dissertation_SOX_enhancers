#!/bin/bash

#$ -N chipQC
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=256G
#$ -pe sharedmem 12
#$ -l h_rt=72:00:00

export PATH=/exports/eddie/scratch/s1749179/anaconda/envs/stats/bin:$PATH

Rscript chipQC.R
