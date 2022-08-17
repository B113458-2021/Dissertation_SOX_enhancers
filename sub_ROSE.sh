#!/bin/bash

#$ -N ROSE
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -pe sharedmem 16
#$ -l h_rt=05:00:00

CONFIG=$1
IDS=$2

source $CONFIG

### ROSE: Rank Ordering of Super-Enhancers
### Based on https://github.com/stjude/ROSE

NARROWPEAK=`head -n $SGE_TASK_ID $IDS | cut -f 1| tail -n 1`
BAMFILE=`head -n $SGE_TASK_ID $IDS | cut -f 2| tail -n 1`
CONTROL=`head -n $SGE_TASK_ID $IDS | cut -f 3| tail -n 1`

PATHTO=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/src/scripts/ROSE
PYTHONPATH=$PATHTO/lib
export PYTHONPATH
export PATH=$PATH:$PATHTO/bin

cd $ROSE_IN

if [ ! -f ${NARROWPEAK%_peaks.narrowPeak}.gff ]; then
    cat $NARROWPEAK | sort -k1,1 -k2,2n | awk -F\\t '{print $1 "\t" NR "\t\t" $2 "\t" $3 "\t\t.\t\t" NR}' > ${NARROWPEAK%_peaks.narrowPeak}.gff
else 
    echo "Processed gff file already exist."
fi

if [ ! -f $ROSE_OUT/${NARROWPEAK%_peaks.narrowPeak}_Plot_points.png ]; then

    if [ $CONTROL != 'NA' ]; then
        ROSE_main.py -g hg38 -i ${NARROWPEAK%_peaks.narrowPeak}.gff -r $BAMFILE -c $CONTROL --custom $ROSE_ANNOTATION -o $ROSE_OUT

    else

        ROSE_main.py -g hg38 -i ${NARROWPEAK%_peaks.narrowPeak}.gff -r $BAMFILE --custom $ROSE_ANNOTATION -o $ROSE_OUT
    fi

fi

