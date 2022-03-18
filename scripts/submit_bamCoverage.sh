#!/bin/bash

## This script uses deeptools suite to generate coverage tracks.
#$ -N bamCov
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=24:00:00
#$ -pe sharedmem 16

IDS=$1
CHIP_BAMS_DIR=$2
COVERAGE_DIR=$3

## Deeptools need python 3.65
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/py365/bin:$PATH

## Take sorted bams and generate bigwig tracks:
SAMPLE_ID=`head -n $SGE_TASK_ID $IDS | cut -f 1 | tail -n 1`

## For paired ends ends
bamCoverage \
   -b $CHIP_BAMS_DIR/${SAMPLE_ID}.sorted.bam \
   -o $COVERAGE_DIR/${SAMPLE_ID}.sorted.bw \
   --outFileFormat bigwig \
   --binSize 10 \
   --scaleFactor 10 \
   --smoothLength 105 \
   --normalizeUsing CPM \
   --effectiveGenomeSize 2913022398 \
   --ignoreForNormalization chrX chrM chrY \
   --extendReads
