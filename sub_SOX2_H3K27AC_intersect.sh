#!/bin/bash

#$ -N SOX2H3K27AC_intersect
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -pe sharedmem 16
#$ -l h_rt=12:00:00

CONFIG=$1
IDS=$2

source $CONFIG
export PATH=$PATH:/exports/eddie/scratch/s1749179/sox_cancer_enhancers/src/scripts/homer/.//bin/

SAMPLE_NAME=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 1`
SOX2=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 2`
H3K27AC=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 3`

cd $MACS2_PEAKS

echo "Intersecting SOX2 and H3K27ac data for ${SAMPLE_NAME}..."

bedtools intersect -a $SOX2 -b $H3K27AC -wa | cut -f 1,2,3 > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps.bed
bedtools intersect -a $H3K27AC -b $SOX2 -wa | cut -f 1,2,3 >> ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps.bed
sort -k1,1 -k2,2n ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps.bed | bedtools merge -i stdin | grep -vsi random | grep -vsi Un | grep -vsi chrM | grep -vsi alt > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged.bed
rm ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps.bed
bedtools intersect -a $SOX2 -b $H3K27AC -wa | cut -f 1,2,3 > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion.bed

echo "Annotating putative enhancer regions..."
annotatePeaks.pl ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged.bed hg38 > ${HOMER_ANNOTATIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged.annotated.txt
annotatePeaks.pl ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion.bed hg38 > ${HOMER_ANNOTATIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion.annotated.txt

echo "Filtering putative enhancers into intergenic, promoter and intragenic..."
cat ${HOMER_ANNOTATIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged.annotated.txt | grep -vsi peakid | grep -si Intergenic | cut -f 2-4 | sort -k1 > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged_intergenic.bed
cat ${HOMER_ANNOTATIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged.annotated.txt | grep -vsi peakid | grep -si TSS | cut -f 2-4 | sort -k1 > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged_promoter.bed
cat ${HOMER_ANNOTATIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged.annotated.txt | grep -vsi peakid | grep -vsi Intergenic | grep -vsi TSS | cut -f 2-4 | sort -k1 > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged_intragenic.bed

cat ${HOMER_ANNOTATIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion.annotated.txt | grep -vsi peakid | grep -si Intergenic | cut -f 2-4 | sort -k1 > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion_intergenic.bed
cat ${HOMER_ANNOTATIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion.annotated.txt | grep -vsi peakid | grep -si TSS | cut -f 2-4 | sort -k1 > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion_promoter.bed
cat ${HOMER_ANNOTATIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion.annotated.txt | grep -vsi peakid | grep -vsi Intergenic | grep -vsi TSS | cut -f 2-4 | sort -k1 > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion_intragenic.bed

echo "Extracting fasta sequences..."
bedtools getfasta -fi $REFERENCE -bed ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged.bed > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged.fasta
bedtools getfasta -fi $REFERENCE -bed ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged_intergenic.bed > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged_intergenic.fasta
bedtools getfasta -fi $REFERENCE -bed ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged_promoter.bed > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged_promoter.fasta
bedtools getfasta -fi $REFERENCE -bed ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged_intragenic.bed > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_merged_intragenic.fasta
bedtools getfasta -fi $REFERENCE -bed ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion.bed > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion.fasta
bedtools getfasta -fi $REFERENCE -bed ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion_intergenic.bed > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion_intergenic.fasta
bedtools getfasta -fi $REFERENCE -bed ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion_promoter.bed > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion_promoter.fasta
bedtools getfasta -fi $REFERENCE -bed ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion_intragenic.bed > ${INTERSECTIONS}/${SAMPLE_NAME}_SOX2_H3K27AC_overlaps_SOXregion_intragenic.fasta
