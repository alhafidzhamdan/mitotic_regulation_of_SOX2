
#!/bin/bash

## This script takes in commands for deeptool's plotHeatmap function

#$ -N all_peaks
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=100G
#$ -l h_rt=12:00:00
    
CONFIG=$1
ANALYSIS_NAME=$2

### GENOMIC_ELEMENT points to a bed path of genomic elements
source $CONFIG

## Deeptools need python 3.65
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/py365/bin:$PATH

BW_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/coverage/merged/normalised
MATRIX_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/downstream/deeptools/heatmap/computeMatrix
HEATMAP_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/downstream/deeptools/heatmap/figures

#### Define dba between phases regions:
ALL_DBA_PEAKS=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/downstream/deeptools/heatmap/deeptools_bed/all_peaks_between_phases.bed

#### Compute matrix:
if [[ ! -f ${MATRIX_DIR}/${ANALYSIS_NAME}.matrix.gz ]];
then
    computeMatrix reference-point \
        --referencePoint center \
        -R $ALL_DBA_PEAKS \
        -S ${BW_DIR}/*.normalised.merged.bw \
        -b 2500 -a 2500 \
        --binSize 10 \
        --missingDataAsZero \
        --blackListFileName $BLACKLISTED_PEAKS_MM10 \
        -p 32 \
        --scale 10 \
        --outFileName ${MATRIX_DIR}/${ANALYSIS_NAME}.matrix.gz \
        --outFileNameMatrix ${MATRIX_DIR}/${ANALYSIS_NAME}.matrix.tab \
        --outFileSortedRegions ${MATRIX_DIR}/${ANALYSIS_NAME}.matrix.bed
fi

#### Plot heatmap:
plotHeatmap \
    -m ${MATRIX_DIR}/${ANALYSIS_NAME}.matrix.gz \
    -out ${HEATMAP_DIR}/${ANALYSIS_NAME}.heatmap.pdf \
    --heatmapHeight 21 \
    --heatmapWidth 3 \
    --legendLocation none \
    --regionsLabel "All peaks" \
    --yAxisLabel "CPM" \
    --xAxisLabel "Peak distance (bp)" \
    --verbose

    
    
    
    
    
    
