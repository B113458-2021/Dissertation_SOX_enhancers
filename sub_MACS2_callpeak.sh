#!/bin/bash

# To run this script, do 
# qsub -t 1-n sub_MACS2_callpeak.sh <CONFIG> <THRESHOLD>


#$ -N MACS2
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=12G
#$ -l h_rt=02:00:00
#$ -pe sharedmem 12


CONFIG=$1
THRESHOLD=$2

source $CONFIG

INPUT=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 5)
BAM_NAME=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 7)
SAMPLE_NAME=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 1)
CHIP_TYPE=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 3)

cd $MACS2_PEAKS

if [ $INPUT == 'YES' ];
then
	if [ $CHIP_TYPE == 'H3K27AC' ];
	then
		macs2 callpeak \
		-t ${BAMS}/${BAM_NAME}.final.bam \
		-c ${BAMS}/${SAMPLE_NAME}_INPUT.final.bam \
		-n ${BAM_NAME}_q${THRESHOLD}_broad \
		-q ${THRESHOLD} \
		-g hs \
		--broad
	else
		macs2 callpeak \
			-t ${BAMS}/${BAM_NAME}.final.bam \
			-c ${BAMS}/${SAMPLE_NAME}_INPUT.final.bam \
			-n ${BAM_NAME}_q${THRESHOLD} \
			-q ${THRESHOLD} \
			-g hs
	fi

else
	if [ $CHIP_TYPE == 'H3K27AC' ];
	then

		macs2 callpeak \
			-t ${BAMS}/${BAM_NAME}.final.bam \
			-n ${BAM_NAME}_q${THRESHOLD}_broad   \
			-q ${THRESHOLD} \
			-g hs \
			--broad
	else
		macs2 callpeak \
			-t ${BAMS}/${BAM_NAME}.final.bam \
			-n ${BAM_NAME}_q${THRESHOLD} \
			-q ${THRESHOLD} \
			-g hs
	fi
fi

