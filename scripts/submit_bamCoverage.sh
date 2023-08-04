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

## For paired ends ends
## Standalone:
## CHIP
bamCoverage \
   -b $CHIP_BAMS_DIR/${SAMPLE_ID}-ChIP.merged.sorted.bam  \
   -o $COVERAGE_DIR/${SAMPLE_ID}-ChIP.merged.rpkm.bw \
   --outFileFormat bigwig \
   --binSize 10 \
   --scaleFactor 10 \
   --smoothLength 105 \
   --normalizeUsing RPKM \
   --effectiveGenomeSize 2730871774 \
   --ignoreForNormalization chrX chrM chrY \
   --extendReads

## INPUT
bamCoverage \
   -b $CHIP_BAMS_DIR/${SAMPLE_ID}-Input.merged.sorted.bam \
   -o $COVERAGE_DIR/${SAMPLE_ID}-Input.merged.rpkm.bw \
   --outFileFormat bigwig \
   --binSize 10 \
   --scaleFactor 10 \
   --smoothLength 105 \
   --normalizeUsing RPKM \
   --effectiveGenomeSize 2730871774 \
   --ignoreForNormalization chrX chrM chrY \
   --extendReads

## Normalised against input file
BLACKLISTED_PEAKS_MM10=/exports/igmm/eddie/Glioblastoma-WGS/resources/mm10-blacklist.v2.bed.gz

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






##for rep in `seq 1 3`
##do
##bamCompare \
##    -b1 $CHIP_BAMS_DIR/${SAMPLE_ID}-ChIP-${rep}.sorted.bam \
##    -b2 $CHIP_BAMS_DIR/${SAMPLE_ID}-Input-${rep}.sorted.bam \
##    -o ${COVERAGE_DIR}/${SAMPLE_ID}_${rep}.normalised.bw \
##    --outFileFormat bigwig \
##    -p 32 \
##    --binSize 10 \
##    --smoothLength 105 \
##    --normalizeUsing CPM \
##    --effectiveGenomeSize 2730871774 \
##    --blackListFileName $BLACKLISTED_PEAKS_MM10 \
##    --ignoreForNormalization chrX chrM chrY \
##    --scaleFactorsMethod None
##done





## Need to specify extension length D for single ends
#####D=`head -n $SGE_TASK_ID $IDS | cut -f 2 | tail -n 1`
##bamCompare \
##    -b1 $CHIP_BAMS_DIR/${SAMPLE_ID}_H3K27AC.final.bam \
##    -b2 $CHIP_BAMS_DIR/${SAMPLE_ID}_INPUT.final.bam \
##    -o ${SAMPLE_ID}_H3K27ac.normalised.bw \
##    --outFileFormat bigwig \
##    --binSize 10 \
##    --smoothLength 105 \
##    --normalizeUsing CPM \
##    --effectiveGenomeSize 2913022398 \
##    --blackListFileName $ENCODE_BLACKLIST \
##    --ignoreForNormalization chrX chrM chrY \
##    --extendReads $D \
##    --scaleFactorsMethod None




## For single ends
##bamCoverage \
##   -b $CHIP_BAMS_DIR/${SAMPLE_ID}_H3K27AC.sorted.bam \
##   -o $COVERAGE_DIR/${SAMPLE_ID}_H3K27AC.sorted.bw \
##   --outFileFormat bigwig \
##   --binSize 10 \
##   --scaleFactor 10 \
##   --smoothLength 105 \
##   --normalizeUsing CPM \
##   --effectiveGenomeSize 2913022398 \
##  --ignoreForNormalization chrX chrM chrY \
##   --extendReads $D

#######




## Archive
## For SOX9
##bamCoverage \
##   -b $CHIP_BAMS_DIR/${SAMPLE_ID}_Sox9_sorted_merged.bam \
##   -o $COVERAGE_DIR/${SAMPLE_ID}_Sox9_sorted_merged.bw \
##   --outFileFormat bigwig \
##   --binSize 10 \
##   --scaleFactor 10 \
##   --smoothLength 105 \
##   --normalizeUsing CPM \
##   --effectiveGenomeSize 2913022398 \
##   --ignoreForNormalization chrX chrM chrY \
##   --extendReads 
## Plot for both Chip and input files:
##CHIP_REP1=$CHIP_BAMS_DIR/${SAMPLE_ID}_${CHIP}_Rep1.final.bam
##CHIP_REP2=$CHIP_BAMS_DIR/${SAMPLE_ID}_${CHIP}_Rep2.final.bam
##INPUT_REP1=$CHIP_BAMS_DIR/${SAMPLE_ID}_Input_Rep1.final.bam
##INPUT_REP2=$CHIP_BAMS_DIR/${SAMPLE_ID}_Input_Rep2.final.bam

##cd $CHIP_BAMS_DIR
## Merge technical replicates:

##if [[ ! -f $CHIP_BAMS_DIR/${SAMPLE_ID}_${CHIP}_merged.final.bam ]]
##then
##samtools merge -f $CHIP_BAMS_DIR/${SAMPLE_ID}_${CHIP}_merged.final.bam $CHIP_REP1 $CHIP_REP2
##samtools merge -f $CHIP_BAMS_DIR/${SAMPLE_ID}_Input_merged.final.bam $INPUT_REP1 $INPUT_REP2
##samtools index -b $CHIP_BAMS_DIR/${SAMPLE_ID}_${CHIP}_merged.final.bam
##samtools index -b $CHIP_BAMS_DIR/${SAMPLE_ID}_Input_merged.final.bam
##fi

## Generate coverage/bedgraph track:
##cd ${CHIP_DEEPTOOLS}/coverage/bedgraph

## For Chip 
##if [[ ! -f ${SAMPLE_ID}_${CHIP}.bedgraph ]]
##then
##echo "Building coverage track for chip file..."
##bamCoverage \
##    -b $CHIP_BAMS_DIR/${SAMPLE_ID}_${CHIP}_merged.final.bam \
##    -o ${SAMPLE_ID}_${CHIP}.bedgraph \
##    --outFileFormat bedgraph \
##    --binSize 10 \
##    --scaleFactor 10 \
##    --smoothLength 105 \
##    --normalizeUsing CPM \
##    --effectiveGenomeSize 2913022398 \
##    --ignoreForNormalization chrX chrM chrY \
##    --extendReads
##else
##echo "Coverage track for CHIP file already generated..."
##fi

## For input
##if [[ ! -f ${SAMPLE_ID}_Input.bedgraph ]]
##then
##echo "Building coverage track for input file..."
##bamCoverage \
##    -b $CHIP_BAMS_DIR/${SAMPLE_ID}_Input_Rep1.final.bam \
##    -o ${SAMPLE_ID}_Input.bedgraph \
##    --outFileFormat bedgraph \
##    --binSize 10 \
##    --scaleFactor 10 \
##    --smoothLength 105 \
##    --normalizeUsing CPM \
##    --effectiveGenomeSize 2913022398 \
##    --ignoreForNormalization chrX chrM chrY \
##    --extendReads
##else
##echo "Coverage track for Input file already generated..."
##fi    


## Optionally prepare merged bam file as above and run pairwise (CHIP vs input) comparison using bamCompare.
##if [[ ! -f ${SAMPLE_ID}_${CHIP}.log2ratio.bedgraph ]]
###then
##echo "Building normalised coverage track for chip file, normalised to input signals..."
##bamCompare \
##    -b1 $CHIP_BAMS_DIR/${SAMPLE_ID}_${CHIP}_merged.final.bam \
##    -b2 $CHIP_BAMS_DIR/${SAMPLE_ID}_Input_Rep1.final.bam \
##    -o ${SAMPLE_ID}_${CHIP}.log2ratio.bedgraph \
##    --outFileFormat bedgraph \
##    --binSize 10 \
##    --smoothLength 105 \
##    --normalizeUsing CPM \
##    --effectiveGenomeSize 2913022398 \
##    --ignoreForNormalization chrX chrM chrY \
##    --extendReads \
##    --scaleFactorsMethod None

##    --scaleFactors 10 \
##else 
##echo "Normalised coverage track already generated..."
###fi

## Subset to certain region:
##if [[ ! -f ${SAMPLE_ID}_${CHIP}.${CHR}.log2ratio.bedgraph ]]
##then
##echo "Subsetting log2ratio file for ${SAMPLE_ID}_${CHIP}.log2ratio.bedgraph file to ${CHR} region..."
##grep -si ${CHR} ${SAMPLE_ID}_${CHIP}.log2ratio.bedgraph > ${SAMPLE_ID}_${CHIP}.${CHR}.log2ratio.bedgraph
##fi

##if [[ ! -f ${SAMPLE_ID}_${CHIP}.${CHR}.bedgraph ]]
##then
####echo "Subsetting bedgraph file for ${SAMPLE_ID}_${CHIP}.bedgraph to ${CHR} region..."
##grep -si ${CHR} ${SAMPLE_ID}_${CHIP}.bedgraph > ${SAMPLE_ID}_${CHIP}.${CHR}.bedgraph
##fi

##if [[ ! -f ${SAMPLE_ID}_Input.${CHR}.bedgraph ]]
##then
##echo "Subsetting input file for ${SAMPLE_ID}_Input to ${CHR}..."
##grep -si ${CHR} ${SAMPLE_ID}_Input.bedgraph > ${SAMPLE_ID}_Input.${CHR}.bedgraph
##fi 

#### For bigwigs:
## Generate coverage/bedgraph track:
##cd ${CHIP_DEEPTOOLS}/coverage/bigwig

## For Chip 
##if [[ ! -f ${SAMPLE_ID}_${CHIP}.bw ]]
##then
##echo "Building coverage track for chip file..."
##bamCoverage \
##    -b $CHIP_BAMS_DIR/${SAMPLE_ID}_${CHIP}_merged.final.bam \
##    -o ${SAMPLE_ID}_${CHIP}.bw \
##    --outFileFormat bigwig \
##    --binSize 10 \
##    --scaleFactor 10 \
##    --smoothLength 105 \
##    --normalizeUsing CPM \
##    --effectiveGenomeSize 2913022398 \
##    --ignoreForNormalization chrX chrM chrY \
##    --extendReads
##else
##echo "Coverage track for CHIP file already generated..."
##fi

## For input
##if [[ ! -f ${SAMPLE_ID}_Input.bw ]]
##then
##echo "Building coverage track for input file..."
##bamCoverage \
##    -b $CHIP_BAMS_DIR/${SAMPLE_ID}_Input_Rep1.final.bam \
##    -o ${SAMPLE_ID}_Input.bw \
##    --outFileFormat bigwig \
##    --binSize 10 \
##    --scaleFactor 10 \
##    --smoothLength 105 \
##    --normalizeUsing CPM \
##    --effectiveGenomeSize 2913022398 \
##    --ignoreForNormalization chrX chrM chrY \
##    --extendReads
##else
##echo "Coverage track for Input file already generated..."
##fi    

##CHIP_BAMS_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/H3K27ac-seq_GBM/reads/Single_end_reads
##CHIP_COVERAGE=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/H3K27ac-seq_GBM/qc/deeptools/coverage/bigwig
##cd $CHIP_COVERAGE

## Optionally prepare merged bam file as above and run pairwise (CHIP vs input) comparison using bamCompare.
##if [[ ! -f ${SAMPLE_ID}_${CHIP}.normalised.bw ]]
##then
##echo "Building normalised coverage track for ${SAMPLE_ID} chip file, normalised to input signals..."
##bamCompare \
#    -b1 $CHIP_BAMS_DIR/${SAMPLE_ID}_H3K27AC.final.bam \
#    -b2 $CHIP_BAMS_DIR/${SAMPLE_ID}_INPUT.final.bam \
#    -o ${SAMPLE_ID}_H3K27ac.normalised.bw \
#    --outFileFormat bigwig \
#    --binSize 10 \
#    --smoothLength 105 \
#    --normalizeUsing CPM \
#    --effectiveGenomeSize 2913022398 \
#    --blackListFileName $ENCODE_BLACKLIST \
#    --ignoreForNormalization chrX chrM chrY \
#    --extendReads 200 \
#    --scaleFactorsMethod None
##else 
##echo "Normalised coverage track already generated..."
#fi

# CHIP_BAMS_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/SOX_TF/bams
# COVERAGE_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/SOX_TF/downstream/coverage/sorted_bams_coverage
# SAMPLE_ID=`head -n $SGE_TASK_ID $IDS | cut -f 1| tail -n 1`
# CHIP=`head -n $SGE_TASK_ID $IDS | cut -f 2| tail -n 1`
# INPUT=`head -n $SGE_TASK_ID $IDS | cut -f 3| tail -n 1`
# D=`head -n $SGE_TASK_ID $IDS | cut -f 4| tail -n 1`
# CHIP_BAM=${CHIP}.final.bam
# INPUT_BAM=${INPUT}.final.bam

## For single-end data need to provide extension value

# bamCoverage \
# -b $CHIP_BAMS_DIR/$CHIP_BAM \
# -o $COVERAGE_DIR/${SAMPLE_ID}.bw \
# --scaleFactor 10 \
# --outFileFormat bigwig \
# --binSize 100 \
# --smoothLength 200 \
# --verbose \
# -p 8 \
# --normalizeUsing RPKM \
# --effectiveGenomeSize 2913022398 \
# --blackListFileName $ENCODE_BLACKLIST \
# --ignoreForNormalization chrX chrM chrY \
# --extendReads $D
# 
# bamCompare \
#    -b1 $CHIP_BAMS_DIR/$CHIP_BAM \
#    -b2 $CHIP_BAMS_DIR/$INPUT_BAM \
#    -o ${CHIP_BAM%_merged.final.bam}.norm.bw \
#    --outFileFormat bigwig \
#    --binSize 100 \
#    --smoothLength 200 \
#    -p 8 \
#    --normalizeUsing CPM \
#    --effectiveGenomeSize 2913022398 \
#    --blackListFileName $ENCODE_BLACKLIST \
#    --ignoreForNormalization chrX chrM chrY \
#    --extendReads 200 \
#    --scaleFactorsMethod None










