#!/bin/bash

# To run this script, do 
# qsub -t 1-n sub_alignment.sh <CONFIG>

#$ -N ALIGNMENT
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=08:00:00
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
        bwa aln -t 16 ${REFERENCE}/hg38_idx ${TRIMMED_FASTQS}/trimmed_${SAMPLE_NAME}_1.fastq.gz > ${BAMS}/${SAMPLE_NAME}_1.sai
        bwa aln -t 16 ${REFERENCE}/hg38_idx ${TRIMMED_FASTQS}/trimmed_${SAMPLE_NAME}_2.fastq.gz > ${BAMS}/${SAMPLE_NAME}_2.sai
        bwa sampe -P ${REFERENCE}/hg38_idx ${BAMS}/${SAMPLE_NAME}_1.sai ${BAMS}/${SAMPLE_NAME}_2.sai ${TRIMMED_FASTQS}/trimmed_${SAMPLE_NAME}_1.fastq.gz ${TRIMMED_FASTQS}/trimmed_${SAMPLE_NAME}_2.fastq.gz | samtools sort -@ 10 -T ${SAMPLE_NAME} -o ${BAMS}/${SAMPLE_NAME}.sorted.bam
        rm ${BAMS}/${SAMPLE_NAME}_1.sai ${BAMS}/${SAMPLE_NAME}_2.sai
else
        bwa aln -t 16 ${REFERENCE}/hg38_idx ${TRIMMED_FASTQS}/trimmed_${SAMPLE_NAME}.fastq.gz > ${BAMS}/${SAMPLE_NAME}.sai
        bwa samse ${REFERENCE}/hg38_idx ${BAMS}/${SAMPLE_NAME}.sai ${TRIMMED_FASTQS}/trimmed_${SAMPLE_NAME}.fastq.gz | samtools sort -@ 10 - -T ${SAMPLE_NAME} -o ${BAMS}/${SAMPLE_NAME}.sorted.bam
        rm ${BAMS}/${SAMPLE_NAME}.sai
fi

samtools index ${BAMS}/${SAMPLE_NAME}.sorted.bam
samtools flagstat ${BAMS}/${SAMPLE_NAME}.sorted.bam > ${BAMS}/${SAMPLE_NAME}.stats.out
