#$ -N bigWigmerge
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=32G
#$ -l h_rt=10:00:00
#$ -pe sharedmem 16


CONFIG=$1

source $CONFIG

bigWigMerge ${COVERAGE}/GSC*H3K27AC.bw ${COVERAGE}/GSC_H3K27AC.merged 
