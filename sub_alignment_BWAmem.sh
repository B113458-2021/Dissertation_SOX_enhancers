#!/bin/bash

# To run this script, do 
# qsub -t 1-n sub_alignment.sh <CONFIG>

#$ -N ALIGNMENT
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=05:00:00
#$ -pe sharedmem 16


CONFIG=$1

source $CONFIG

TYPE=$(tail -n +2 ${SAMPLES} | head -n $SGE_TASK_ID | tail -n 1 | cut -f 4)
SAMPLE_NAME=$(tail -n +2 ${SAMPLES} | head -n $SGE_TASK_ID | tail -n 1 | cut -f 6)

# First need to create a BWA hg38 genome index by running
# bwa index -p ${REFERENCE}/hg38_idx -a bwtsw ${REFERENCE}/hg38.fa.gz


#then align fastq reads to reference using BWA aln
if [ "$TYPE" == "PAIRED" ];
then
        bwa mem -M -t 20 ${REFERENCE}/hg38_idx ${TRIMMED_FASTQS}/trimmed_${SAMPLE_NAME}_1.fastq.gz ${TRIMMED_FASTQS}/trimmed_${SAMPLE_NAME}_2.fastq.gz | samtools sort -@ 10 -T ${SAMPLE_NAME} -o ${BAMS}/${SAMPLE_NAME}.sorted.bam
else
        bwa mem -M -t 20 ${REFERENCE}/hg38_idx ${TRIMMED_FASTQS}/trimmed_${SAMPLE_NAME}.fastq.gz | samtools sort -@ 10 - -T ${SAMPLE_NAME} -o ${BAMS}/${SAMPLE_NAME}.sorted.bam
fi

samtools index ${BAMS}/${SAMPLE_NAME}.sorted.bam
samtools flagstat ${BAMS}/${SAMPLE_NAME}.sorted.bam > ${BAMS}/${SAMPLE_NAME}.stats.out
