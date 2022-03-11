#!/bin/bash

# To run this script, do 
# qsub -t 1-n submit_plotFingerprint.sh CONFIG IDS CHIP BATCH
#
#$ -N plotFingerprint
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=6G
#$ -l h_rt=72:00:00
#$ -pe sharedmem 6 

## QC on filtered bam files using deeptools
## Based on https://deeptools.readthedocs.io/en/develop/index.html
## Generate various metrics useful for quality control.

CONFIG=$1
IDS=$2
###CHIP=$3

SAMPLE_ID=`head -n $SGE_TASK_ID $IDS | tail -n 1`

source $CONFIG

## Deeptools need python 3.65
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/py365/bin:$PATH

CHIP=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/H3K27ac-seq_GBM/bams/${SAMPLE_ID}_H3K27AC.final.bam
INPUT=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/H3K27ac-seq_GBM/bams/${SAMPLE_ID}_INPUT.final.bam

cd /exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/H3K27ac-seq_GBM/qc/deeptools

echo "Plotting fingerprint for ${SAMPLE_ID}..."
plotFingerprint \
    --bamfiles $CHIP $INPUT \
    --plotFile ${SAMPLE_ID}_H3K27ac.fingerprint.png \
    --minMappingQuality 30 --skipZeros \
    --smartLabels --extendReads 150 \
    --region 1 --numberOfSamples 50000 
    






#CHIP_REP1=$CHIP_BAMS_DIR/${SAMPLE_ID}_${CHIP}_Rep1.final.bam
#INPUT_REP1=$CHIP_BAMS_DIR/${SAMPLE_ID}_Input_Rep1.final.bam
#CHIP_REP2=$CHIP_BAMS_DIR/${SAMPLE_ID}_${CHIP}_Rep2.final.bam
##INPUT_REP2=$CHIP_BAMS_DIR/${SAMPLE_ID}_Input_Rep2.final.bam

#cd $CHIP_DEEPTOOLS
#plotFingerprint \
#    --bamfiles $CHIP_REP1 $CHIP_REP2 $INPUT_REP1 \
#    --plotFile ${SAMPLE_ID}_${CHIP}.fingerprint.png \
#    --minMappingQuality 30 --skipZeros \
#    --smartLabels \
#    --outRawCounts ${SAMPLE_ID}_${CHIP}.fingerprint.counts \
#    --region 1 --numberOfSamples 50000 
    
