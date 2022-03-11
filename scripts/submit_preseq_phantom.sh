#!/bin/bash

# To run this script, do 
# qsub -t 1-n submit_preseq_phantom.sh CONFIG IDS
#
#
#$ -N prsqphntm
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=20:00:00
#$ -pe sharedmem 8

CONFIG=$1
IDS=$2

source $CONFIG

SAMPLE_ID=`head -n $SGE_TASK_ID $IDS | tail -n 1`
PRESEQ_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/H3K27ac-seq_GBM/qc/preseq
BAM_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/H3K27ac-seq_GBM/bams
BAM_FILE=${SAMPLE_ID}.sorted.bam

## Phantompeakqual and preseq 
## Need R 3.6, spp package installed.
## This checks for library complexity (need for further sequencing/ probability of missing peaks), and presence of phantom peaks
## Must use unprocessed/unfiltered sorted indexed bams
## Based on https://scilifelab.github.io/courses/ngsintro/1604/labs/chipseq

## Use preseq package to calculate library complexity
preseq c_curve -v -B -o $PRESEQ_DIR/${SAMPLE_ID}_preseq_estimates.txt $BAM_DIR/$BAM_FILE 2> $PRESEQ_DIR/${SAMPLE_ID}_verbose.txt

## NSC and RSC correlates
cd $PHANTOMPEAKQUAL

Rscript run_spp.R \
    -c=$BAM_DIR/$BAM_FILE \
    -savp=$PRESEQ_DIR/${SAMPLE_ID}_xcor.pdf \
    -out=$PRESEQ_DIR/${SAMPLE_ID}_xcor_metrics.txt
    

