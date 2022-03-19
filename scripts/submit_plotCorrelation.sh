#!/bin/bash

# To run this script, do 
# qsub submit_plotCorrelation.sh CONFIG BATCH
#
#$ -N plotCorr
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=120G
#$ -l h_rt=12:00:00

## QC on original sorted bam files using deeptools
## Based on https://deeptools.readthedocs.io/en/develop/index.html
## Generate multiBamSummary for all bam files:

BATCH=$1
CHIP_BAMS_DIR=$2
CHIP_DEEPTOOLS_DIR=$3


## Deeptools need python 3.65
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/py365/bin:$PATH


## Gather all sorted bam files and generate summary of chr1 in a numpy array
cd $CHIP_BAMS_DIR
BAM_FILES=`ls *.final.bam | grep -v / | tr '\n' " "`

if [[ ! -f $CHIP_DEEPTOOLS_DIR/${BATCH}.npz ]]
then
echo "Preparing multibam summary bins..."
multiBamSummary bins \
    --bamfiles $BAM_FILES \
    -o $CHIP_DEEPTOOLS_DIR/${BATCH}.npz
fi

## Plot correlations 
echo "Plotting correlation heatmap..."
plotCorrelation -in $CHIP_DEEPTOOLS_DIR/${BATCH}.npz \
    --removeOutliers \
    --whatToPlot heatmap \
    --skipZeros \
    --corMethod spearman \
    -o $CHIP_DEEPTOOLS_DIR/${BATCH}.corr.spearman.png \
    --outFileCorMatrix $CHIP_DEEPTOOLS_DIR/${BATCH}.corr.spearman.readCounts.tab

## Plot PCA
echo "Plotting PCA..."
plotPCA -in $CHIP_DEEPTOOLS_DIR/${BATCH}.npz \
    --transpose \
    --plotHeight 20 \
    --plotWidth 30 \
     --outFileNameData $CHIP_DEEPTOOLS_DIR/${BATCH}.PCA.transposed.tab \
    -o $CHIP_DEEPTOOLS_DIR/${BATCH}.PCA.png 



echo "Done!"

