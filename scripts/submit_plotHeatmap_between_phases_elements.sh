#!/bin/bash

## This script takes in commands for deeptool's plotHeatmap function, focusing on genomic element of interest ie promoter/introns

#$ -N element_heatmap
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=120G
#$ -l h_rt=10:00:00
    
CONFIG=$1
REGION=$2
LABEL=$3
ANALYSIS_NAME=$4
VARIANT=$5

### GENOMIC_ELEMENT points to a bed path of genomic elements
source $CONFIG

## Deeptools need python 3.65
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/py365/bin:$PATH

#### Plot heatmaps of diffbind bed for phase differences for each variant (exclude mitotic enriched as very few regions)
BW_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/coverage/merged/normalised
MATRIX_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/downstream/deeptools/heatmap/computeMatrix
HEATMAP_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/downstream/deeptools/heatmap/figures

if [[ ! -f ${MATRIX_DIR}/${VARIANT}_${ANALYSIS_NAME}.matrix.gz ]];
then
    computeMatrix reference-point \
        --referencePoint center \
        -R $REGION \
        -S ${BW_DIR}/${VARIANT}-Interphase.normalised.merged.bw ${BW_DIR}/${VARIANT}-Mitotic.normalised.merged.bw \
        -b 5000 -a 5000 \
        --binSize 10 \
        --missingDataAsZero \
        --blackListFileName $BLACKLISTED_PEAKS_MM10 \
        -p 32 \
        --scale 10 \
        --samplesLabel ${VARIANT}_Interphase ${VARIANT}_Mitosis \
        --outFileName ${MATRIX_DIR}/${VARIANT}_${ANALYSIS_NAME}.matrix.gz \
        --outFileNameMatrix ${MATRIX_DIR}/${VARIANT}_${ANALYSIS_NAME}.matrix.tab \
        --outFileSortedRegions ${MATRIX_DIR}/${VARIANT}_${ANALYSIS_NAME}.matrix.bed
fi

plotHeatmap \
    -m ${MATRIX_DIR}/${VARIANT}_${ANALYSIS_NAME}.matrix.gz \
    -out ${HEATMAP_DIR}/${VARIANT}_${ANALYSIS_NAME}.heatmap.pdf \
    --heatmapHeight 21 \
    --heatmapWidth 3 \
    --legendLocation none \
    --colorMap PuBu \
    --regionsLabel $LABEL \
    --yAxisLabel "CPM" \
    --xAxisLabel "Peak distance (bp)" \
    --verbose
    
    
    
    
    
