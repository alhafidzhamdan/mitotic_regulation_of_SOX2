#!/bin/bash

## This script takes in commands for deeptool's plotHeatmap function centred around CDS

#$ -N cds_heatmap
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=150G
#$ -l h_rt=4:00:00
    
CONFIG=$1
ANALYSIS_NAME=$2
VARIANT=$3

source $CONFIG

## Deeptools need python 3.65
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/py365/bin:$PATH

#### Plot heatmaps of diffbind bed for phase differences for each variant (exclude mitotic enriched as very few regions)
BW_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/coverage/merged/normalised
MATRIX_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/downstream/deeptools/heatmap/computeMatrix
HEATMAP_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/downstream/deeptools/heatmap/figures

#### CDS:
CDS_GTF=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/references/gencode.vM25.basic.annotation.gtf.gz

if [[ ! -f ${MATRIX_DIR}/${VARIANT}_${ANALYSIS_NAME}.matrix.gz ]];
then
    computeMatrix scale-regions \
        -R $CDS_GTF \
        --metagene \
        --exonID gene_id \
        -S ${BW_DIR}/${VARIANT}-Interphase.normalised.merged.bw ${BW_DIR}/${VARIANT}-Mitotic.normalised.merged.bw \
        -b 2500 -a 2500 \
        --binSize 10 \
        --regionBodyLength 3000 \
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
    --colorMap PuOr \
    --regionsLabel "Genes" \
    --yAxisLabel "CPM" \
    --xAxisLabel "Peak distance (bp)" \
    --verbose

    
    
    
    
    
    
