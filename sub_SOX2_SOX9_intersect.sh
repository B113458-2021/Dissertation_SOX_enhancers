#!/bin/bash

#$ -N SOX2_SOX9_intersect
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=01:00:00


CONFIG=$1
IDS=$2


source $CONFIG
export PATH=$PATH:/exports/eddie/scratch/s1749179/sox_cancer_enhancers/src/scripts/homer/.//bin/

SAMPLE_NAME=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 1`
SOX2=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 2`
SOX9=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 3`

cd $MACS2_PEAKS

echo "Intersecting SOX2 and SOX9 data for ${SAMPLE_NAME}..."

## Find SOX2 and SOX9 co-binding regions:
bedtools intersect -a $SOX2 -b $SOX9 | grep -vsi random | grep -vsi Un | grep -vsi chrM | grep -vsi alt | grep -vsi fix | sortBed -i - | cut -f 1,2,3 > ${INTERSECTING_REGIONS}/${SAMPLE_NAME}_Sox2_Sox9_overlap.bed
bedtools getfasta -fi ${REFERENCE}/hg38.fa -bed ${INTERSECTING_REGIONS}/${SAMPLE_NAME}_Sox2_Sox9_overlap.bed > ${INTERSECTING_FASTAS}/${SAMPLE_NAME}_Sox2_Sox9_overlap.fasta

## Find SOX2 region only
bedtools intersect -a $SOX2 -b $SOX9 -v | grep -vsi random | grep -vsi Un | grep -vsi chrM | grep -vsi alt | grep -vsi fix | sortBed -i - | cut -f 1,2,3 > ${INTERSECTING_REGIONS}/${SAMPLE_NAME}_Sox2_only.bed
bedtools getfasta -fi ${REFERENCE}/hg38.fa -bed ${INTERSECTING_REGIONS}/${SAMPLE_NAME}_Sox2_only.bed > ${INTERSECTING_FASTAS}/${SAMPLE_NAME}_Sox2_only.fasta

## Find SOX9 region only
bedtools intersect -a $SOX9 -b $SOX2 -v | grep -vsi random | grep -vsi Un | grep -vsi chrM | grep -vsi alt | grep -vsi fix | sortBed -i - | cut -f 1,2,3 > ${INTERSECTING_REGIONS}/${SAMPLE_NAME}_Sox9_only.bed
bedtools getfasta -fi ${REFERENCE}/hg38.fa -bed ${INTERSECTING_REGIONS}/${SAMPLE_NAME}_Sox9_only.bed > ${INTERSECTING_FASTAS}/${SAMPLE_NAME}_Sox9_only.fasta
