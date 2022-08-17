#!/bin/bash

# To run this script, do 
# qsub -t 1-n sub_alignment.sh <CONFIG>

#$ -N BAMQC
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=05:00:00
#$ -pe sharedmem 8

CONFIG=$1

source $CONFIG

cd $BAMS

SAMPLE_NAME=$(tail -n +2 ${SAMPLES} | head -n $SGE_TASK_ID | tail -n 1 | cut -f 6)

## Remove multi-mapped and supplementary reads and retain only uniquely mapped reads.
## Code based on https://bioinformatics.stackexchange.com/questions/508/obtaining-uniquely-mapped-reads-from-bwa-mem-alignment

samtools view -h -F 4 ${SAMPLE_NAME}.sorted.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' > ${SAMPLE_NAME}.unique.sorted.bam
samtools sort -@ 10 -n -o ${SAMPLE_NAME}.sorted.unique.sorted.bam ${SAMPLE_NAME}.unique.sorted.bam
rm ${SAMPLE_NAME}.unique.sorted.bam
samtools fixmate -m ${SAMPLE_NAME}.sorted.unique.sorted.bam ${SAMPLE_NAME}.fixmate.sorted.unique.sorted.bam
rm ${SAMPLE_NAME}.sorted.unique.sorted.bam
samtools sort -@ 10 -o ${SAMPLE_NAME}.position.fixmate.sorted.unique.sorted.bam ${SAMPLE_NAME}.fixmate.sorted.unique.sorted.bam
rm ${SAMPLE_NAME}.fixmate.sorted.unique.sorted.bam
samtools markdup ${SAMPLE_NAME}.position.fixmate.sorted.unique.sorted.bam ${SAMPLE_NAME}.markdup.position.fixmate.sorted.unique.sorted.bam
rm ${SAMPLE_NAME}.position.fixmate.sorted.unique.sorted.bam 

## Filter blacklisted regions

bedtools intersect -v -a ${SAMPLE_NAME}.markdup.position.fixmate.sorted.unique.sorted.bam -b $BLACKLISTED_PEAKS > ${SAMPLE_NAME}.filtered.markdup.position.fixmate.sorted.unique.sorted.bam
rm ${SAMPLE_NAME}.markdup.position.fixmate.sorted.unique.sorted.bam
samtools sort -@ 10 -o ${SAMPLE_NAME}.final.bam ${SAMPLE_NAME}.filtered.markdup.position.fixmate.sorted.unique.sorted.bam
rm ${SAMPLE_NAME}.filtered.markdup.position.fixmate.sorted.unique.sorted.bam
samtools index ${SAMPLE_NAME}.final.bam
