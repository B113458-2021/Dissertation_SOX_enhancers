#!/bin/bash


# To run this script, do 
# qsub -t 1-n sub_TrimFilt.sh <CONFIG>

#$ -N TRIMFILTER
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=05:00:00
#$ -pe sharedmem 8

CONFIG=$1

source $CONFIG

TYPE=$(tail -n +2 ${SAMPLES} | head -n $SGE_TASK_ID | tail -n 1 | cut -f 4)
SAMPLE_NAME=$(tail -n +2 ${SAMPLES} | head -n $SGE_TASK_ID | tail -n 1 | cut -f 6)

if [ "$TYPE" == "PAIRED" ];
then
	#trim and filter the raw reads to remove adapters and discard low quality reads
	fastp -i ${FASTQS}/${SAMPLE_NAME}_1.fastq.gz -I ${FASTQS}/${SAMPLE_NAME}_2.fastq.gz -o ${TRIMMED_FASTQS}/trimmed_${SAMPLE_NAME}_1.fastq.gz -O ${TRIMMED_FASTQS}/trimmed_${SAMPLE_NAME}_2.fastq.gz -V --trim_poly_g --detect_adapter_for_pe -h ${READ_QC}/${SAMPLE_NAME}.fastp.html -j ${READ_QC}/${SAMPLE_NAME}.fastp.json
	#run quality control again on the trimmed data
	fastqc ${TRIMMED_FASTQS}/trimmed_${SAMPLE_NAME}_1.fastq.gz ${TRIMMED_FASTQS}/trimmed_${SAMPLE_NAME}_2.fastq.gz -t 2 -o $READ_QC
else
	fastp -i ${FASTQS}/${SAMPLE_NAME}.fastq.gz -o ${TRIMMED_FASTQS}/trimmed_${SAMPLE_NAME}.fastq.gz -V --trim_poly_g -h ${READ_QC}/${SAMPLE_NAME}.fastp.html -j ${READ_QC}/${SAMPLE_NAME}.fastp.json

	fastqc ${TRIMMED_FASTQS}/trimmed_${SAMPLE_NAME}.fastq.gz -o $READ_QC
fi

