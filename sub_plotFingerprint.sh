#!/bin/bash

# To run this script, do 
# qsub -t 1-n submit_plotFingerprint.sh CONFIG
#
#$ -N plotFingerprint
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=05:00:00
#$ -pe sharedmem 8 

## QC on filtered bam files using deeptools
## Based on https://deeptools.readthedocs.io/en/develop/index.html
## Generate various metrics useful for quality control.

CONFIG=$1

source $CONFIG

SAMPLE_NAME=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 1)
CHIP_TYPE=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 3)
TYPE=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 4)
INPUT=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 5)
REPLICATES=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 6)
BAM_NAME=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 7)
D=$(tail -n +2 $BAM_SAMPLES | head -n $SGE_TASK_ID | tail -n 1 | cut -f 8)

cd $FINGERPRINT_OUT

if [ $INPUT == 'YES' ];
then

    if [ $REPLICATES == '1' ];
    then
        if [ $TYPE == 'SINGLE' ];
        then
            echo "Plotting fingerprint for ${BAM_NAME}..."
            plotFingerprint \
                --bamfiles ${BAMS}/${BAM_NAME}.final.bam ${BAMS}/${SAMPLE_NAME}_INPUT.final.bam \
                --plotFile ${BAM_NAME}.fingerprint.png \
                --minMappingQuality 30 --skipZeros \
                --smartLabels --extendReads $D \
                --outRawCounts ${BAM_NAME}.fingerprint.counts \
                --region 1 --numberOfSamples 50000

        elif [ $TYPE == 'PAIRED' ];
        then
            echo "Plotting fingerprint for ${BAM_NAME}..."
            plotFingerprint \
                --bamfiles ${BAMS}/${BAM_NAME}.final.bam ${BAMS}/${SAMPLE_NAME}_INPUT.final.bam \
                --plotFile ${BAM_NAME}.fingerprint.png \
                --minMappingQuality 30 --skipZeros \
                --smartLabels \
                --outRawCounts ${BAM_NAME}.fingerprint.counts \
                --region 1 --numberOfSamples 50000
        fi

    elif [ $REPLICATES == '2' ];
    then
        if [ $TYPE == 'SINGLE' ];
        then
            echo "Plotting fingerprint for ${BAM_NAME}..."
            plotFingerprint \
                --bamfiles ${BAMS}/${SAMPLE_NAME}_R1_${CHIP_TYPE}.final.bam ${BAMS}/${SAMPLE_NAME}_R2_${CHIP_TYPE}.final.bam ${BAMS}/${SAMPLE_NAME}_INPUT.final.bam \
                --plotFile ${BAM_NAME}.fingerprint.png \
                --minMappingQuality 30 --skipZeros \
                --smartLabels --extendReads $D \
                --outRawCounts ${BAM_NAME}.fingerprint.counts \
                --region 1 --numberOfSamples 50000

        elif [ $TYPE == 'PAIRED' ];
        then
            echo "Plotting fingerprint for ${BAM_NAME}..."
            plotFingerprint \
                --bamfiles ${BAMS}/${SAMPLE_NAME}_R1_${CHIP_TYPE}.final.bam ${BAMS}/${SAMPLE_NAME}_R2_${CHIP_TYPE}.final.bam ${BAMS}/${SAMPLE_NAME}_INPUT.final.bam \
                --plotFile ${BAM_NAME}.fingerprint.png \
                --minMappingQuality 30 --skipZeros \
                --smartLabels \
                --outRawCounts ${BAM_NAME}.fingerprint.counts \
                --region 1 --numberOfSamples 50000
        fi
    fi
else
    echo "No input control for ${BAM_NAME}."
    if [ $REPLICATES == '1' ];
    then
        if [ $TYPE == 'SINGLE' ];
        then
            echo "Plotting fingerprint for ${BAM_NAME}..."
            plotFingerprint \
                --bamfiles ${BAMS}/${BAM_NAME}.final.bam \
                --plotFile ${BAM_NAME}.fingerprint.png \
                --minMappingQuality 30 --skipZeros \
                --smartLabels --extendReads $D \
                --outRawCounts ${BAM_NAME}.fingerprint.counts \
                --region 1 --numberOfSamples 50000

        elif [ $TYPE == 'PAIRED' ];
        then
            echo "Plotting fingerprint for ${BAM_NAME}..."
            plotFingerprint \
                --bamfiles ${BAMS}/${BAM_NAME}.final.bam \
                --plotFile ${BAM_NAME}.fingerprint.png \
                --minMappingQuality 30 --skipZeros \
                --smartLabels \
                --outRawCounts ${BAM_NAME}.fingerprint.counts \
                --region 1 --numberOfSamples 50000
        fi

    elif [ $REPLICATES == '2' ];
    then
        if [ $TYPE == 'SINGLE' ];
        then
            echo "Plotting fingerprint for ${BAM_NAME}..."
            plotFingerprint \
                --bamfiles ${BAMS}/${SAMPLE_NAME}_R1_${CHIP_TYPE}.final.bam ${BAMS}/${SAMPLE_NAME}_R2_${CHIP_TYPE}.final.bam \
                --plotFile ${BAM_NAME}.fingerprint.png \
                --minMappingQuality 30 --skipZeros \
                --smartLabels --extendReads $D \
                --outRawCounts ${BAM_NAME}.fingerprint.counts \
                --region 1 --numberOfSamples 50000

        elif [ $TYPE == 'PAIRED' ];
        then
            echo "Plotting fingerprint for ${BAM_NAME}..."
            plotFingerprint \
                --bamfiles ${BAMS}/${SAMPLE_NAME}_R1_${CHIP_TYPE}.final.bam ${BAMS}/${SAMPLE_NAME}_R2_${CHIP_TYPE}.final.bam \
                --plotFile ${BAM_NAME}.fingerprint.png \
                --minMappingQuality 30 --skipZeros \
                --smartLabels \
                --outRawCounts ${BAM_NAME}.fingerprint.counts \
                --region 1 --numberOfSamples 50000
        fi
    fi

fi
