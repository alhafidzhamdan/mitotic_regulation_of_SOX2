#!/bin/bash

#$ -N plotHeatmapProfile
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=100:00:00
####$ -pe sharedmem 8 
    
CONFIG=$1
SAMPLE_ID=$2
source $CONFIG

## Deeptools need python 3.65
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/py365/bin:$PATH

BED=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/SOX_TF/downstream/beds/SOX_SE_for_deeptools.bed
##PROMOTER=${SAMPLE_ID}_Sox2_Sox9_merged_overlap_base.annotated.promoter.bed
##INTERGENIC=${SAMPLE_ID}_Sox2_Sox9_merged_overlap_base.annotated.intergenic.bed

CHIP_DEEPTOOLS=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/SOX_TF/qc/deeptools
OUTPUT_DIR=${CHIP_DEEPTOOLS}/computeMatrix/${SAMPLE_ID}_SOX_SE_targets

if [ ! -d $OUTPUT_DIR ]; then
    mkdir -p $OUTPUT_DIR && cd $OUTPUT_DIR
else
    cd $OUTPUT_DIR
fi

##if [ ! -f ${SAMPLE_ID}_co_bound_SE_targets.matrix.gz ]; then
computeMatrix reference-point \
    --referencePoint center \
    -R $BED \
    -S ${CHIP_DEEPTOOLS}/coverage/bigwig/*_Sox2.normalised.bw \
    -b 1000 -a 1000 \
    --binSize 10 \
    --skipZeros --missingDataAsZero \
    --blackListFileName $ENCODE_BLACKLIST \
    -p 6 \
    -o ${SAMPLE_ID}_co_bound_SE_targets.Sox2.matrix.gz \
    --outFileNameMatrix ${SAMPLE_ID}_co_bound_SE_targets.Sox2.matrix.tab \
    --outFileSortedRegions ${SAMPLE_ID}_co_bound_SE_targets.Sox2.matrix.bed \
    --verbose 
##fi

plotHeatmap \
    -m ${SAMPLE_ID}_co_bound_SE_targets.Sox2.matrix.gz \
    -out ${SAMPLE_ID}_co_bound_SE_targets.Sox2.heatmap.png \
    --heatmapHeight 15 \
    --legendLocation none \
    --verbose 

computeMatrix reference-point \
    --referencePoint center \
    -R $BED \
    -S ${CHIP_DEEPTOOLS}/coverage/bigwig/*_Sox9.normalised.bw \
    -b 1000 -a 1000 \
    --binSize 10 \
    --skipZeros --missingDataAsZero \
    --blackListFileName $ENCODE_BLACKLIST \
    -p 6 \
    -o ${SAMPLE_ID}_co_bound_SE_targets.Sox9.matrix.gz \
    --outFileNameMatrix ${SAMPLE_ID}_co_bound_SE_targets.Sox9.matrix.tab \
    --outFileSortedRegions ${SAMPLE_ID}_co_bound_SE_targets.Sox9.matrix.bed \
    --verbose 
##fi

plotHeatmap \
    -m ${SAMPLE_ID}_co_bound_SE_targets.Sox9.matrix.gz \
    -out ${SAMPLE_ID}_co_bound_SE_targets.Sox9.heatmap.png \
    --heatmapHeight 15 \
    --legendLocation none \
    --verbose 
    
    
#plotProfile -m ${CHIP_SEQ_ID}.matrix.gz \
#    -out ${CHIP_SEQ_ID}.pergroup.Sox2_centre_250.profile.png \
#    --perGroup 
    
    
