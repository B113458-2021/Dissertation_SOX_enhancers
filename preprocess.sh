#!/bin/bash 

FASTQ_PATH=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/raw_data/GSC_H3K27ac

tail -n +2 ${FASTQ_PATH}/GSC_datasets_H3K27AC.txt | while read line || [ -n "$line" ];
do
	echo $line
	OLDNAME=$(echo "$line" | awk 'BEGIN {FS="\t"} {print $5}')
	NEWNAME=$(echo "$line" | awk 'BEGIN {FS="\t"} {print $6}')
	TYPE=$(echo "$line" | awk 'BEGIN {FS="\t"} {print $4}')

	echo $OLDNAME
	echo $NEWNAME
	echo $TYPE
	
	if [ "$TYPE" == "PAIRED" ]; 
	then
		mv ${FASTQ_PATH}/${OLDNAME}_R1.fastq.gz $FASTQ_PATH/${NEWNAME}_1.fastq.gz
		mv ${FASTQ_PATH}/${OLDNAME}_R2.fastq.gz $FASTQ_PATH/${NEWNAME}_2.fastq.gz
	else
		mv ${FASTQ_PATH}/${OLDNAME}.fastq.gz $FASTQ_PATH/${NEWNAME}.fastq.gz
	fi
done
