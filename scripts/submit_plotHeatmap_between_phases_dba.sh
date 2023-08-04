#!/bin/bash

## This script takes in commands for deeptool's plotHeatmap function
## Example run (in directory "logs"): qsub ../../../scripts/fame_gbm/mitotic_bookmarking/submit_plotHeatmap_between_phases_dba.sh ../../../scripts/fame_gbm/config.sh between_phases_all WT

#$ -N btwn_phases_heatmap
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=150G
#$ -l h_rt=4:00:00
    
CONFIG=$1
ANALYSIS_NAME=$2
VARIANT=$3

### GENOMIC_ELEMENT points to a bed path of genomic elements
source $CONFIG

## Deeptools need python 3.65
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/py365/bin:$PATH

#### Plot heatmaps of diffbind bed for phase differences for each variant (exclude mitotic enriched as very few regions)
BW_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/coverage/merged/normalised
MATRIX_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/downstream/deeptools/heatmap/computeMatrix
HEATMAP_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/downstream/deeptools/heatmap/figures

#### Define dba between phases regions:
BETWEEN_PHASES_DBA=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/diffbind/extended_2q1vsq5/dba_beds/between_phases
INTERPHASE_ENRICHED=${BETWEEN_PHASES_DBA}/${VARIANT}_interphase_enriched_FC0.bed
SHARED=${BETWEEN_PHASES_DBA}/${VARIANT}_shared_FC0.bed
MITOSIS_ENRICHED=${BETWEEN_PHASES_DBA}/${VARIANT}_mitosis_enriched_FC0.bed

#### Compute matrix:
if [[ ! -f ${MATRIX_DIR}/${VARIANT}_${ANALYSIS_NAME}.matrix.gz ]];
then
    computeMatrix reference-point \
        --referencePoint center \
        -R $INTERPHASE_ENRICHED $SHARED \
        -S ${BW_DIR}/${VARIANT}-Interphase.normalised.merged.bw ${BW_DIR}/${VARIANT}-Mitotic.normalised.merged.bw \
        -b 2500 -a 2500 \
        --binSize 100 \
        --missingDataAsZero \
        --blackListFileName $BLACKLISTED_PEAKS_MM10 \
        -p 32 \
        --scale 10 \
        --samplesLabel ${VARIANT}_Interphase ${VARIANT}_Mitosis \
        --outFileName ${MATRIX_DIR}/${VARIANT}_${ANALYSIS_NAME}.matrix.gz \
        --outFileNameMatrix ${MATRIX_DIR}/${VARIANT}_${ANALYSIS_NAME}.matrix.tab \
        --outFileSortedRegions ${MATRIX_DIR}/${VARIANT}_${ANALYSIS_NAME}.matrix.bed
fi

#### Plot heatmap:
plotHeatmap \
    -m ${MATRIX_DIR}/${VARIANT}_${ANALYSIS_NAME}.matrix.gz \
    -out ${HEATMAP_DIR}/${VARIANT}_${ANALYSIS_NAME}.heatmap.pdf \
    --heatmapHeight 21 \
    --heatmapWidth 3 \
    --legendLocation none \
    --regionsLabel "Interphase enriched" "Shared" \
    --yAxisLabel "CPM" \
    --xAxisLabel "Peak distance (bp)" \
    --verbose

    
    
    
    
    
    
