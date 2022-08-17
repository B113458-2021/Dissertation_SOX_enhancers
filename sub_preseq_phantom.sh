#!/bin/bash

# To run this script, do 
# qsub -t 1-n submit_preseq_phantom.sh CONFIG
#
#
#$ -N prsqphntm
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=03:00:00
#$ -pe sharedmem 8

CONFIG=$1

source $CONFIG

SAMPLE_NAME=$(tail -n +2 ${SAMPLES} | head -n $SGE_TASK_ID | tail -n 1 | cut -f 6)
BAM_FILE=${SAMPLE_NAME}.sorted.bam

## Phantompeakqual and preseq 
## Need R 3.6, spp package installed.
## This checks for library complexity (need for further sequencing/ probability of missing peaks), and presence of phantom peaks
## Must use unprocessed/unfiltered sorted indexed bams
## Based on https://scilifelab.github.io/courses/ngsintro/1604/labs/chipseq

## Use preseq package to calculate library complexity
preseq c_curve -v -B -o $PRESEQ_OUT/${SAMPLE_NAME}_preseq_estimates.txt $BAMS/$BAM_FILE 2> $PRESEQ_OUT/${SAMPLE_NAME}_verbose.txt

## NSC and RSC correlates
cd $PHANTOMPEAKQUAL

Rscript run_spp.R \
    -c=$BAMS/$BAM_FILE \
    -savp=$PRESEQ_OUT/${SAMPLE_NAME}_xcor.pdf \
    -out=$PRESEQ_OUT/${SAMPLE_NAME}_xcor_metrics.txt
