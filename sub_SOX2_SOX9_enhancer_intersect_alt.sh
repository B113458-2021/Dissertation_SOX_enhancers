#!/bin/bash

#$ -N SOX2_SOX9_enhancer_intersect
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=1:00:00


CONFIG=$1

source $CONFIG

#### Generate consensus SOX2/SOX9 and SOX2 enhancers for GSCs
echo "Generating consensus SOX2/SOX9 and SOX2 enhancers for GSCs..."

### Sox2/Sox9
bedtools multiinter -i ${INTERSECTING_REGIONS}/E17_Sox2_Sox9_overlap.bed ${INTERSECTING_REGIONS}/E21_Sox2_Sox9_overlap.bed ${INTERSECTING_REGIONS}/E27_Sox2_Sox9_overlap.bed ${INTERSECTING_REGIONS}/E28_Sox2_Sox9_overlap.bed ${INTERSECTING_REGIONS}/E31_Sox2_Sox9_overlap.bed ${INTERSECTING_REGIONS}/E34_Sox2_Sox9_overlap.bed ${INTERSECTING_REGIONS}/E37_Sox2_Sox9_overlap.bed > ${INTERSECTING_REGIONS}/GSC_SOX2_SOX9_Peaks_counts.bed 

## Extract only those peaks that are present/overlapped by at least 5 (majority) of the GSC cell lines --> GSC SOX2/SOX9 consensus peaks
awk '$4 >= 5 {print}' ${INTERSECTING_REGIONS}/GSC_SOX2_SOX9_Peaks_counts.bed | sortBed -i - | awk '{print $1"\t"$2"\t"$3}' > ${INTERSECTING_REGIONS}/GSC_SOX2_SOX9_consensusUnmerged.bed
bedtools merge -d 100 -i ${INTERSECTING_REGIONS}/GSC_SOX2_SOX9_consensusUnmerged.bed | awk '($3-$2>=50) {print}'| uniq > ${INTERSECTING_REGIONS}/GSC_SOX2_SOX9_consensusPeaks.bed

### Sox2 only
bedtools multiinter -i ${MACS2_PEAKS}/E17_SOX2_q0.05_peaks.narrowPeak ${MACS2_PEAKS}/E21_SOX2_q0.05_peaks.narrowPeak ${MACS2_PEAKS}/E27_SOX2_q0.05_peaks.narrowPeak ${MACS2_PEAKS}/E28_SOX2_q0.05_peaks.narrowPeak ${MACS2_PEAKS}/E31_SOX2_q0.05_peaks.narrowPeak ${MACS2_PEAKS}/E34_SOX2_q0.05_peaks.narrowPeak ${MACS2_PEAKS}/E37_SOX2_q0.05_peaks.narrowPeak > ${INTERSECTING_REGIONS}/GSC_SOX2_Peaks_counts.bed

## Extract only those peaks that are present/overlapped by at least 5 (majority) of the GSC cell lines --> GSC SOX2 consensus peaks
awk '$4 >= 5 {print}' ${INTERSECTING_REGIONS}/GSC_SOX2_Peaks_counts.bed | grep -vsi random | grep -vsi Un | grep -vsi chrM | grep -vsi alt | grep -vsi fix | sortBed -i - | awk '{print $1"\t"$2"\t"$3}' > ${INTERSECTING_REGIONS}/GSC_SOX2_consensusUnmerged.bed
bedtools merge -d 100 -i ${INTERSECTING_REGIONS}/GSC_SOX2_consensusUnmerged.bed | awk '($3-$2>=50) {print}'| uniq > ${INTERSECTING_REGIONS}/GSC_SOX2_consensusPeaks.bed

#### Generate consensus SOX2 enhancers for SQUAMOUS
echo "Generating consensus SOX2 enhancers for Squamous cancer cell lines..."

bedtools multiinter -i ${MACS2_PEAKS}/KYSE140_SOX2_q0.05_peaks.narrowPeak ${MACS2_PEAKS}/LK2_SOX2_q0.05_peaks.narrowPeak ${MACS2_PEAKS}/NCI-H520_SOX2_q0.05_peaks.narrowPeak ${MACS2_PEAKS}/TE5_SOX2_q0.05_peaks.narrowPeak ${MACS2_PEAKS}/SCC25_SOX2_q0.05_peaks.narrowPeak > ${INTERSECTING_REGIONS}/Squamous_SOX2_Peaks_counts.bed

## Extract only those peaks that are present/overlapped by at least 3 (majority) of the Squmaous cell lines --> Squamous SOX2 consensus peaks

awk '$4 >= 3 {print}' ${INTERSECTING_REGIONS}/Squamous_SOX2_Peaks_counts.bed | grep -vsi random | grep -vsi Un | grep -vsi chrM | grep -vsi alt | grep -vsi fix | sortBed -i - | awk '{print $1"\t"$2"\t"$3}' > ${INTERSECTING_REGIONS}/Squamous_SOX2_consensusUnmerged.bed
bedtools merge -d 100 -i ${INTERSECTING_REGIONS}/Squamous_SOX2_consensusUnmerged.bed | awk '($3-$2>=50) {print}'| uniq > ${INTERSECTING_REGIONS}/Squamous_SOX2_consensusPeaks.bed


#### GSC vs. SQUAMOUS comparison
echo "Comparing GSC enhancers with squamous enhancers..."
## Extract GSC SOX2/SOX9 consensus peaks that are present/overlap the squamous cancer cell lines SOX2 peaks --> GSC/Squamous overlapping SOX2/SOX9 peaks

bedtools intersect -wa \
	-a ${INTERSECTING_REGIONS}/GSC_SOX2_SOX9_consensusPeaks.bed \
	-b ${INTERSECTING_REGIONS}/Squamous_SOX2_consensusPeaks.bed | uniq | sortBed -i - > ${INTERSECTING_REGIONS}/GSC_squamous_SOX2_SOX9_overlapping_final.bed

## Extract GSC SOX2/SOX9 consensus peaks that are completely unique to GSCs --> GSC unique SOX2/SOX9 peaks

bedtools intersect -v \
	-a ${INTERSECTING_REGIONS}/GSC_SOX2_SOX9_consensusPeaks.bed \
	-b ${INTERSECTING_REGIONS}/Squamous_SOX2_consensusPeaks.bed | uniq | sortBed -i - > ${INTERSECTING_REGIONS}/GSC_unique_SOX2_SOX9_final.bed

## Extract Squamous SOX2 consensus peaks that are completely unique to Squamous (when compared to GSC SOX2/SOX9) --> Squamous unique SOX2 peaks (SOX2/SOX9)

bedtools intersect -v \
	-a ${INTERSECTING_REGIONS}/Squamous_SOX2_consensusPeaks.bed \
	-b  ${INTERSECTING_REGIONS}/GSC_SOX2_SOX9_consensusPeaks.bed | uniq | sortBed -i - > ${INTERSECTING_REGIONS}/Squamous_unique_SOX2_SOX9_final.bed

## Extract SOX2 consensus peaks that are present/overlap between the GSC cell lines and squamous cancer cell lines -> GSC/Squamous overlapping SOX2 peaks

bedtools intersect \
	-a ${INTERSECTING_REGIONS}/GSC_SOX2_consensusPeaks.bed \
	-b ${INTERSECTING_REGIONS}/Squamous_SOX2_consensusPeaks.bed | uniq | sortBed -i - > ${INTERSECTING_REGIONS}/GSC_squamous_SOX2_overlapping_final.bed

## Extract GSC SOX2 consensus peaks that are completely unique to GSCs --> GSC unique SOX2 peaks

bedtools intersect -v \
	-a ${INTERSECTING_REGIONS}/GSC_SOX2_consensusPeaks.bed \
	-b ${INTERSECTING_REGIONS}/Squamous_SOX2_consensusPeaks.bed | uniq | sortBed -i - > ${INTERSECTING_REGIONS}/GSC_unique_SOX2_final.bed

## Extract Squamous SOX2 consensus peaks that are completely unique to Squamous --> Squamous unique SOX2 peaks

bedtools intersect -v \
	-a ${INTERSECTING_REGIONS}/Squamous_SOX2_consensusPeaks.bed \
	-b ${INTERSECTING_REGIONS}/GSC_SOX2_consensusPeaks.bed | uniq | sortBed -i - > ${INTERSECTING_REGIONS}/Squamous_unique_SOX2_final.bed

## Extract the fasta sequences for all comparisons
echo "Extracting fasta sequences..."
bedtools getfasta -fi ${REFERENCE}/hg38.fa -bed ${INTERSECTING_REGIONS}/GSC_squamous_SOX2_SOX9_overlapping_final.bed > ${INTERSECTING_FASTAS}/GSC_squamous_SOX2_SOX9_overlapping_final.fasta
bedtools getfasta -fi ${REFERENCE}/hg38.fa -bed ${INTERSECTING_REGIONS}/GSC_unique_SOX2_SOX9_final.bed > ${INTERSECTING_FASTAS}/GSC_unique_SOX2_SOX9_final.fasta
bedtools getfasta -fi ${REFERENCE}/hg38.fa -bed ${INTERSECTING_REGIONS}/Squamous_unique_SOX2_SOX9_final.bed > ${INTERSECTING_FASTAS}/Squamous_unique_SOX2_SOX9_final.fasta
bedtools getfasta -fi ${REFERENCE}/hg38.fa -bed ${INTERSECTING_REGIONS}/GSC_squamous_SOX2_overlapping_final.bed > ${INTERSECTING_FASTAS}/GSC_squamous_SOX2_overlapping_final.fasta
bedtools getfasta -fi ${REFERENCE}/hg38.fa -bed ${INTERSECTING_REGIONS}/GSC_unique_SOX2_final.bed > ${INTERSECTING_FASTAS}/GSC_unique_SOX2_final.fasta
bedtools getfasta -fi ${REFERENCE}/hg38.fa -bed ${INTERSECTING_REGIONS}/Squamous_unique_SOX2_final.bed > ${INTERSECTING_FASTAS}/Squamous_unique_SOX2_final.fasta
