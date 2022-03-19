#!/bin/bash

# To run this script, do 
# qsub -t 1-n submit_filterChip-seq.sh CONFIG IDS
#
# CONFIG is the path to the file scripts/config.sh which contains environment variables set to
# commonly used paths and files in the script
# IDS is a list of sample ids, one per line, where tumor and normal samples are designated by
# the addition of a T or an N to the sample id.
# TYPE is either sample name or input
#
#$ -N filterChip-seq
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=24:00:00
#$ -pe sharedmem 8

CONFIG=$1
IDS=$2
BAM_DIR=$3

source $CONFIG

SAMPLE_ID=`head -n $SGE_TASK_ID $IDS | tail -n 1`

cd $BAM_DIR

################################################################################################
### BWA actually follows the SAM spec and reports Phred scores as MAPQ values.  
### The calculation is based on the number of optimal (best) alignments found, 
### as well as the number of sub-optimal alignments combined with the 
### Phred scores of the bases which differ between the optimal and sub-optimal alignments.
### https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/

REP=${SAMPLE_ID}.sorted.bam
UNIQUE_REP=${SAMPLE_ID}.unique.sorted.bam
SORTED_UNIQUE_REP=${SAMPLE_ID}.sorted.unique.sorted.bam
FIXMATE_REP=${SAMPLE_ID}.fixmate.sorted.unique.sorted.bam
POSITION_FIXMATE_REP=${SAMPLE_ID}.position.fixmate.sorted.unique.sorted.bam
MARKDUP_REP=${SAMPLE_ID}.markdup.position.fixmate.sorted.unique.sorted.bam
FILTERED_REP=${SAMPLE_ID}.filtered.markdup.position.fixmate.sorted.unique.sorted.bam
FINAL_REP=${SAMPLE_ID}.final.bam

## Remove multi-mapped and supplementary reads and retain only uniquely mapped reads.
## Based on https://bioinformatics.stackexchange.com/questions/508/obtaining-uniquely-mapped-reads-from-bwa-mem-alignment
samtools view -h -F 4 $REP | grep -v -e 'XA:Z:' -e 'SA:Z:' > $UNIQUE_REP
samtools sort -@ 10 -n -o $SORTED_UNIQUE_REP $UNIQUE_REP
samtools fixmate -m $SORTED_UNIQUE_REP $FIXMATE_REP
samtools sort -@ 10 -o $POSITION_FIXMATE_REP $FIXMATE_REP
samtools markdup $POSITION_FIXMATE_REP $MARKDUP_REP

## Filter blacklisted regions
## Peaks downloaded from https://github.com/Boyle-Lab/Blacklist/tree/master/lists
bedtools intersect -v -a $MARKDUP_REP -b $BLACKLISTED_PEAKS_MM10 > $FILTERED_REP
samtools sort -@ 10 -o $FINAL_REP $FILTERED_REP
samtools index $FINAL_REP
