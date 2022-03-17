#!/bin/bash

# To run this script, do 
# qsub -t 1-n submit_bwa.sh CONFIG IDS
#
# CONFIG is the path to the file scripts/config.sh which contains environment variables set to
# commonly used paths and files in the script
# IDS is a list of sample ids, one per line, where tumor and normal samples are designated by
# the addition of a T or an N to the sample id.
# TYPE is either sample name or input
#
#$ -N BWA
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=16G
#$ -l h_rt=90:00:00
#$ -pe sharedmem 16

CONFIG=$1
IDS=$2
READ_DIR=$3
BAM_DIR=$4

source $CONFIG

SAMPLE_ID=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 8`

## Index ref genome already downloaded.
## Align fastq reads 

##if paired-end reads
cd $READ_DIR
echo "Aligning fastq file for ${SAMPLE_ID} ... "
bwa mem -M -t 20 $BWA_REF ${SAMPLE_ID}_R1_001.fastq.gz ${SAMPLE_ID}_R2_001.fastq.gz | samtools sort -@ 10 - -T ${SAMPLE_ID} -o $BAM_DIR/${SAMPLE_ID}.sorted.bam
samtools index $BAM_DIR/${SAMPLE_ID}.sorted.bam
samtools flagstat $BAM_DIR/${SAMPLE_ID}.sorted.bam > $BAM_DIR/${SAMPLE_ID}.stats.out
