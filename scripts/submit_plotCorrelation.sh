#!/bin/bash

# To run this script, do 
# qsub submit_plotCorrelation.sh CONFIG BATCH
#
#$ -N plotCorr
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=32G
#$ -l h_rt=72:00:00

## QC on original sorted bam files using deeptools
## Based on https://deeptools.readthedocs.io/en/develop/index.html
## Generate multiBamSummary for all bam files:

CONFIG=$1
BATCH=$2

source $CONFIG

## Deeptools need python 3.65
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/py365/bin:$PATH

## Gather all sorted bam files and generate summary in a numpy array
cd $CHIP_BAMS_DIR
BAM_FILES=`ls *chr1.final.bam | grep -v / | tr '\n' " "`

if [[ ! -f $CHIP_DEEPTOOLS/${BATCH}.npz ]]
then
echo "Preparing multibam summary bins..."
multiBamSummary bins \
    --bamfiles $BAM_FILES \
    -o $CHIP_DEEPTOOLS/${BATCH}.npz
fi

## Plot correlations 
echo "Plotting correlation heatmap..."
plotCorrelation -in $CHIP_DEEPTOOLS/${BATCH}.npz \
    --removeOutliers \
    --whatToPlot heatmap \
    --skipZeros \
    --corMethod spearman \
    -o $CHIP_DEEPTOOLS/${BATCH}.corr.spearman.png \
    --outFileCorMatrix $CHIP_DEEPTOOLS/${BATCH}.corr.spearman.readCounts.tab

## Plot PCA
echo "Plotting PCA..."
plotPCA -in $CHIP_DEEPTOOLS/${BATCH}.npz \
    --transpose \
    --plotHeight 20 \
    --plotWidth 30 \
     --outFileNameData $CHIP_DEEPTOOLS/${BATCH}.PCA.transposed.tab \
    -o $CHIP_DEEPTOOLS/${BATCH}.PCA.png 


