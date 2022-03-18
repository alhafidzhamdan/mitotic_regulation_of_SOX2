#!/bin/bash

# This script detects and filters out adaptor sequences from sequencing data
# https://github.com/OpenGene/fastp
#
#$ -N fastp
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -l h_rt=10:00:00

export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/snakemake/bin:$PATH

IDS=$1
RAW_OUTPUT_DIR=$2
CLEANED_OUTPUT_DIR=$3
QC_DIR=$4
SUFFIX=$5

SAMPLE_ID=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 1`

echo "Automatically detecting adaptor sequences and removing them..."
fastp -i ${RAW_OUTPUT_DIR}/${SAMPLE_ID}_R1${SUFFIX}.fastq.gz -I ${RAW_OUTPUT_DIR}/${SAMPLE_ID}_R2${SUFFIX}.fastq.gz \
     -o ${CLEANED_OUTPUT_DIR}/${SAMPLE_ID}_R1.fastq.gz -O ${CLEANED_OUTPUT_DIR}/${SAMPLE_ID}_R2.fastq.gz \
     -h ${QC_DIR}/${SAMPLE_ID}.fastp.html

echo "Performing fastqc on cleaned fastqs..."
module load igmm/apps/FastQC/0.11.9
fastqc -t 20 ${CLEANED_OUTPUT_DIR}/${SAMPLE_ID}_R1.fastq.gz
fastqc -t 20 ${CLEANED_OUTPUT_DIR}/${SAMPLE_ID}_R2.fastq.gz


