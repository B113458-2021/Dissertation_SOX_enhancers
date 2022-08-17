#!/bin/bash

# To run this script, do 
# qsub -t 1-n sub_IDR.sh <CONFIG> <IDS>

#$ -N IDR
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=01:30:00
#$ -pe sharedmem 8

CONFIG=$1
IDS=$2

source $CONFIG

## Use IDR tool to measure quality of replicates.
## Based on https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/07_handling-replicates-idr.html


SAMPLE_NAME=$(head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 1)
CHIP_TYPE=$(head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 2)



cd $IDR


## Run MACS2 with less stringent cut off for each replicate

################################################################################################################
## Rep 1:

REP1=${SAMPLE_NAME}_R1_${CHIP_TYPE}.final.bam
INPUT=${SAMPLE_NAME}_R1_INPUT.final.bam

macs2 callpeak \
    -t ${BAMS}/${REP1} \
    -c ${BAMS}/${INPUT} \
    -n ${SAMPLE_NAME}_R1_${CHIP_TYPE} \
    -f BAM \
    -B \
    -p 1e-3

## Sort peak by -log10(p-value)

sort -k8,8nr ${SAMPLE_NAME}_R1_${CHIP_TYPE}_peaks.narrowPeak > ${SAMPLE_NAME}_R1_${CHIP_TYPE}_peaks.sorted.narrowPeak

################################################################################################################
## Rep 2: 

REP2=${SAMPLE_NAME}_R2_${CHIP_TYPE}.final.bam

macs2 callpeak \
	-t ${BAMS}/${REP2} \
	-c ${BAMS}/${INPUT} \
	-n ${SAMPLE_NAME}_R2_${CHIP_TYPE} \
	-f BAM \
	-B \
	-p 1e-3

sort -k8,8nr ${SAMPLE_NAME}_R2_${CHIP_TYPE}_peaks.narrowPeak > ${SAMPLE_NAME}_R2_${CHIP_TYPE}_peaks.sorted.narrowPeak


################################################################################################################
## IDR on both runs using default Rank method (in this case signal.value)

idr \
    --samples ${SAMPLE_NAME}_R1_${CHIP_TYPE}_peaks.sorted.narrowPeak ${SAMPLE_NAME}_R2_${CHIP_TYPE}_peaks.sorted.narrowPeak \
    --input-file-type narrowPeak \
    --output-file-type narrowPeak \
    --output-file $IDR/${SAMPLE_NAME}_${CHIP_TYPE}.macs2idr.narrowPeak \
    --plot 


################################################################################################################
