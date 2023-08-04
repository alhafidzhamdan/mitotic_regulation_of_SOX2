#!/bin/bash

## This script takes in commands for deeptool's plotHeatmap function
## Example run (in directory "logs"): qsub ../../../scripts/fame_gbm/mitotic_bookmarking/submit_plotHeatmap_between_variants_dba.sh ../../../scripts/fame_gbm/config.sh between_variants_all WT

#$ -N btwn_variants_heatmap
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=150G
#$ -l h_rt=4:00:00
    
CONFIG=$1
ANALYSIS_NAME=$2
PHASE=$3
DBA_PHASE=$4
VARIANT_A=$5
VARIANT_B=$6

### GENOMIC_ELEMENT points to a bed path of genomic elements
source $CONFIG

## Deeptools need python 3.65
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/py365/bin:$PATH

#### Plot heatmaps of diffbind bed for phase differences for each variant (exclude mitotic enriched as very few regions)
BW_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/coverage/merged/normalised
MATRIX_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/downstream/deeptools/heatmap/computeMatrix
HEATMAP_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/downstream/deeptools/heatmap/figures

#### Define dba regions between variants:
BETWEEN_VARIANTS_DBA=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/diffbind/extended_2q1vsq5/dba_beds/between_variants
A_ENRICHED=${BETWEEN_VARIANTS_DBA}/${DBA_PHASE}_${VARIANT_A}vs${VARIANT_B}_${VARIANT_A}_enriched_FC0.bed
SHARED=${BETWEEN_VARIANTS_DBA}/${DBA_PHASE}_${VARIANT_A}vs${VARIANT_B}_shared_FC0.bed
B_ENRICHED=${BETWEEN_VARIANTS_DBA}/${DBA_PHASE}_${VARIANT_A}vs${VARIANT_B}_${VARIANT_B}_enriched_FC0.bed

#### Compute matrix:
if [[ ! -f ${MATRIX_DIR}/${DBA_PHASE}_${VARIANT_A}vs${VARIANT_B}_${ANALYSIS_NAME}.matrix.gz ]];
then
    computeMatrix reference-point \
        --referencePoint center \
        -R $A_ENRICHED $SHARED $B_ENRICHED \
        -S ${BW_DIR}/${VARIANT_B}-${PHASE}.normalised.merged.bw ${BW_DIR}/${VARIANT_A}-${PHASE}.normalised.merged.bw \
        -b 2500 -a 2500 \
        --binSize 100 \
        --missingDataAsZero \
        --blackListFileName $BLACKLISTED_PEAKS_MM10 \
        -p 32 \
        --scale 10 \
        --samplesLabel ${VARIANT_B}_${DBA_PHASE} ${VARIANT_A}_${DBA_PHASE} \
        --outFileName ${MATRIX_DIR}/${DBA_PHASE}_${VARIANT_A}vs${VARIANT_B}_${ANALYSIS_NAME}.matrix.gz \
        --outFileNameMatrix ${MATRIX_DIR}/${DBA_PHASE}_${VARIANT_A}vs${VARIANT_B}_${ANALYSIS_NAME}.matrix.tab \
        --outFileSortedRegions ${MATRIX_DIR}/${DBA_PHASE}_${VARIANT_A}vs${VARIANT_B}_${ANALYSIS_NAME}.matrix.bed
fi

#### Plot heatmap:
plotHeatmap \
    -m ${MATRIX_DIR}/${DBA_PHASE}_${VARIANT_A}vs${VARIANT_B}_${ANALYSIS_NAME}.matrix.gz \
    -out ${HEATMAP_DIR}/${DBA_PHASE}_${VARIANT_A}vs${VARIANT_B}_${ANALYSIS_NAME}.heatmap.pdf \
    --heatmapHeight 21 \
    --heatmapWidth 3 \
    --legendLocation none \
    --yAxisLabel "CPM" \
    --xAxisLabel "Peak distance (bp)" \
    --verbose

    
    
    
    
    
    
