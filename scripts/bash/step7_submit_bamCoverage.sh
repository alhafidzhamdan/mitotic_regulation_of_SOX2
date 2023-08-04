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

## mmg10 size 2730871774

## Normalised against input file
BLACKLISTED_PEAKS_MM10=/exports/igmm/eddie/Glioblastoma-WGS/resources/mm10-blacklist.v2.bed.gz

## Both technical replicates of TF and input bam files were merged before running this script

bamCompare \
    -b1 $CHIP_BAMS_DIR/${SAMPLE_ID}-ChIP.merged.sorted.bam \
    -b2 $CHIP_BAMS_DIR/${SAMPLE_ID}-Input.merged.sorted.bam \
    -o ${COVERAGE_DIR}/${SAMPLE_ID}.normalised.merged.rpkm.bw \
    --outFileFormat bigwig \
    -p 32 \
    --binSize 10 \
    --smoothLength 105 \
    --normalizeUsing RPKM \
    --effectiveGenomeSize 2730871774 \
    --blackListFileName $BLACKLISTED_PEAKS_MM10 \
    --ignoreForNormalization chrX chrM chrY \
    --scaleFactorsMethod None
