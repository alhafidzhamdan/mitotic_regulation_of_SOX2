#!/bin/bash

# To run this script, do 
# qsub -t 1-n submit_macs2IDR.sh CONFIG IDS CHIP BATCH
# IDS points to sample root name ie G328 (only)
#
#$ -N macs2IDR
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=42:00:00
#$ -pe sharedmem 8

CONFIG=$1
IDS=$2
CHIP=$3
BATCH=$4

source $CONFIG

unset MODULEPATH
. /etc/profile.d/modules.sh


module load igmm/apps/MACS2/2.1.1  

SAMPLE_ID=`head -n $SGE_TASK_ID $IDS | tail -n 1`

## Use IDR tool to measure quality of replicates.
## Based on https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/07_handling-replicates-idr.html

cd $CHIP_IDR/macs_tmp

## Run MACS2 with less stringent cut off for each replicate

################################################################################################################
## Rep 1:

CHIP_REP1=${SAMPLE_ID}_${CHIP}_Rep1.final.bam
INPUT_REP1=${SAMPLE_ID}_Input_Rep1.final.bam

if [[ ! -f ${SAMPLE_ID}_${CHIP}_Rep1_peaks.sorted.narrowPeak && ${SAMPLE_ID}_${CHIP}_Rep2_peaks.sorted.narrowPeak ]]
then
macs2 callpeak \
    -t $CHIP_BAMS_DIR/$CHIP_REP1 \
    -c $CHIP_BAMS_DIR/$INPUT_REP1 \
    -n ${SAMPLE_ID}_${CHIP}_Rep1 \
    -f BAM \
    -B \
    -p 1e-3
    
## Sort peak by -log10(p-value)
sort -k8,8nr ${SAMPLE_ID}_${CHIP}_Rep1_peaks.narrowPeak > ${SAMPLE_ID}_${CHIP}_Rep1_peaks.sorted.narrowPeak

################################################################################################################
## Rep 2: 

CHIP_REP2=${SAMPLE_ID}_${CHIP}_Rep2.final.bam
INPUT_REP1=${SAMPLE_ID}_Input_Rep1.final.bam

macs2 callpeak \
    -t $CHIP_BAMS_DIR/$CHIP_REP2 \
    -c $CHIP_BAMS_DIR/$INPUT_REP1 \
    -n ${SAMPLE_ID}_${CHIP}_Rep2 \
    -f BAM \
    -B \
   -p 1e-3
    
## Sort peak by -log10(p-value)
sort -k8,8nr ${SAMPLE_ID}_${CHIP}_Rep2_peaks.narrowPeak > ${SAMPLE_ID}_${CHIP}_Rep2_peaks.sorted.narrowPeak

fi

################################################################################################################
## IDR on both runs using default Rank method (in this case signal.value)

idr \
    --samples ${SAMPLE_ID}_${CHIP}_Rep1_peaks.sorted.narrowPeak ${SAMPLE_ID}_${CHIP}_Rep2_peaks.sorted.narrowPeak \
    --input-file-type narrowPeak \
    --output-file-type narrowPeak \
    --output-file $CHIP_IDR/${SAMPLE_ID}_${CHIP}.macs2idr.narrowPeak \
    --plot 


################################################################################################################

