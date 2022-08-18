#!/usr/bin/bash

#configuration file containing direcotry and file locations for all scripts

##paths:
export PATH=/exports/eddie/scratch/s1749179/anaconda/envs/chipseq/bin:$PATH

##directories:
SAMPLES=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/raw_data/GSC_datasets.txt
FASTQS=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/raw_data
READ_QC=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/results/QC/FASTQC
TRIMMED_FASTQS=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/processed_data/trimmed_reads
BAMS=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/processed_data/bams
REFERENCE=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/resources/reference
BLACKLISTED_PEAKS=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/resources/encode4_GRCh38_blacklist.bed.gz
MACS2_PEAKS=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/results/MACS2_peaks
IDR=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/results/QC/IDR
BAM_SAMPLES=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/resources/BAM_samples.txt
PRESEQ_OUT=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/results/QC/preseq
FINGERPRINT_OUT=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/results/QC/fingerprint
PHANTOMPEAKQUAL=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/src/scripts/phantompeakqualtools
COVERAGE=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/results/coverage_tracks
ROSE_IN=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/results/ROSE/INPUT
ROSE_OUT=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/results/ROSE/OUTPUT
ROSE_ANNOTATION=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/src/scripts/ROSE/annotation/hg38_refseq.ucsc
INTERSECTING_REGIONS=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/processed_data/intersections/regions
INTERSECTING_FASTAS=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/processed_data/intersections/fastas
HOMER_ANNOTATIONS=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/results/homer_annotations
OVERLAP_STATS=/exports/eddie/scratch/s1749179/sox_cancer_enhancers/results/overlaps_stats

