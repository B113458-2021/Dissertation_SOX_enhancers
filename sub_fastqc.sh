#!/bin/bash

# To run this script, do 
# qsub -t 1-n sub_fastqc.sh <CONFIG>

#$ -N FASTQC
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=00:30:00
#$ -pe sharedmem 8

CONFIG=$1

source $CONFIG

TYPE=$(tail -n +2 ${SAMPLES} | head -n $SGE_TASK_ID | tail -n 1 | cut -f 4)
SAMPLE_NAME=$(tail -n +2 ${SAMPLES} | head -n $SGE_TASK_ID | tail -n 1 | cut -f 6)

if [ "$TYPE" == "PAIRED" ];
then
        fastqc ${FASTQS}/${SAMPLE_NAME}_1.fastq.gz ${FASTQS}/${SAMPLE_NAME}_2.fastq.gz -t 2 -o $READ_QC
else
        fastqc ${FASTQS}/${SAMPLE_NAME}.fastq.gz -o $READ_QC
fi
