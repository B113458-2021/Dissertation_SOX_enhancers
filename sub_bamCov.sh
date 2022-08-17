#!/bin/bash

## This script uses deeptools suite to generate coverage tracks.
#$ -N bamCov
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=10G
#$ -l h_rt=10:00:00
#$ -pe sharedmem 10

CONFIG=$1

source $CONFIG

SAMPLE_NAME=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 1)
TYPE=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 4)
INPUT=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 5)
BAM_NAME=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 7)
D=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 8)




## For paired ends ends
## Standalone:
## CHIP

if [ $TYPE == 'PAIRED' ];
then
	bamCoverage \
	   -b ${BAMS}/${BAM_NAME}.final.bam  \
	   -o ${COVERAGE}/${BAM_NAME}.bw \
	   --outFileFormat bigwig \
	   --binSize 10 \
	   --scaleFactor 10 \
	   --smoothLength 105 \
	   --normalizeUsing CPM \
	   --effectiveGenomeSize 3049315783 \
	   --ignoreForNormalization chrX chrM chrY \
	   --extendReads \
	   --verbose


	if [ $INPUT == 'YES' ];
	then

	   	if [ ! -f "${BAMS}/${SAMPLE_NAME}_INPUT.bw" ];
	   	then
			bamCoverage \
				-b ${BAMS}/${SAMPLE_NAME}_INPUT.final.bam  \
				-o ${COVERAGE}/${SAMPLE_NAME}_INPUT.bw \
				--outFileFormat bigwig \
				--binSize 10 \
				--scaleFactor 10 \
				--smoothLength 105 \
				--normalizeUsing CPM \
				--effectiveGenomeSize 3049315783 \
				--ignoreForNormalization chrX chrM chrY \
				--extendReads \
				--verbose
		fi

			bamCompare \
		    	-b1 ${BAMS}/${BAM_NAME}.final.bam \
		    	-b2 ${BAMS}/${SAMPLE_NAME}_INPUT.final.bam \
		    	-o ${COVERAGE}/${BAM_NAME}.normalised.bw \
		    	--outFileFormat bigwig \
		    	-p 32 \
		    	--binSize 10 \
		    	--smoothLength 105 \
			    --normalizeUsing CPM \
			    --effectiveGenomeSize 3049315783 \
			    --blackListFileName $BLACKLISTED_PEAKS \
			    --ignoreForNormalization chrX chrM chrY \
			    --scaleFactorsMethod None \
			    --verbose
	fi

elif [ $TYPE == 'SINGLE' ];
then

	bamCoverage \
	   -b ${BAMS}/${BAM_NAME}.final.bam  \
	   -o ${COVERAGE}/${BAM_NAME}.bw \
	   --outFileFormat bigwig \
	   --binSize 10 \
	   --scaleFactor 10 \
	   --smoothLength 105 \
	   --normalizeUsing CPM \
	   --effectiveGenomeSize 3049315783 \
	   --ignoreForNormalization chrX chrM chrY \
	   --extendReads $D \
	   --verbose

	
	if [ $INPUT == 'YES' ];
	then

	   	if [ ! -f "${BAMS}/${SAMPLE_NAME}_INPUT.bw" ];
	   	then
			bamCoverage \
				-b ${BAMS}/${SAMPLE_NAME}_INPUT.final.bam  \
				-o ${COVERAGE}/${SAMPLE_NAME}_INPUT.bw \
				--outFileFormat bigwig \
				--binSize 10 \
				--scaleFactor 10 \
				--smoothLength 105 \
				--normalizeUsing CPM \
				--effectiveGenomeSize 3049315783 \
				--ignoreForNormalization chrX chrM chrY \
				--extendReads $D \
				--verbose   
		fi


		bamCompare \
		    	-b1 ${BAMS}/${BAM_NAME}.final.bam \
		    	-b2 ${BAMS}/${SAMPLE_NAME}_INPUT.final.bam \
		    	-o ${COVERAGE}/${BAM_NAME}.normalised.bw \
		    	--outFileFormat bigwig \
		    	-p 32 \
		    	--binSize 10 \
		    	--smoothLength 105 \
			    --normalizeUsing CPM \
			    --effectiveGenomeSize 3049315783 \
			    --blackListFileName $BLACKLISTED_PEAKS \
			    --ignoreForNormalization chrX chrM chrY \
			    --extendReads $D \
			    --scaleFactorsMethod None \
			    --verbose

	fi


fi

