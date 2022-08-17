#!/bin/bash

#$ -N mergeBams
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=02:00:00

CONFIG=$1
IDS=$2
source $CONFIG

echo "Merging and indexing bam files..."

cd $BAMS

SAMPLE_NAME=$(head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 1)
CHIP_TYPE=$(head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 2)

REP1=${SAMPLE_NAME}_R1_${CHIP_TYPE}.final.bam
REP2=${SAMPLE_NAME}_R2_${CHIP_TYPE}.final.bam

samtools merge ${SAMPLE_NAME}_${CHIP_TYPE}.final.bam $REP1 $REP2; samtools index ${SAMPLE_NAME}_${CHIP_TYPE}.final.bam


if [ ! -f "${SAMPLE_NAME}_INPUT.final.bam" ];
then

	#if [ -f "${SAMPLE_NAME}_R2_INPUT.final.bam" ];
	#then
		#INPUT_REP1="${SAMPLE_NAME}_R1_INPUT.final.bam" 
		#INPUT_REP2="${SAMPLE_NAME}_R2_INPUT.final.bam"
		#samtools merge ${SAMPLE_NAME}_INPUT.final.bam $INPUT_REP1 $INPUT_REP2; samtools index ${SAMPLE_NAME}_INPUT.final.bam
	#else
		cp ${SAMPLE_NAME}_R1_INPUT.final.bam ${SAMPLE_NAME}_INPUT.final.bam
		cp ${SAMPLE_NAME}_R1_INPUT.final.bam.bai ${SAMPLE_NAME}_INPUT.final.bam.bai
	#fi
fi
