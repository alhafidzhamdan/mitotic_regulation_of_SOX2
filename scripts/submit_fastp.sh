#!/bin/bash

# This script detects and filters out adaptor sequences from sequencing data
# https://github.com/OpenGene/fastp
# To run this script, do 
# qsub -t 1-N submit_fastp.sh IDS
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
BATCH=$2

PATIENT_ID=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 1`
READ_OUTPUT=/exports/igmm/eddie/Glioblastoma-WGS/RNA-seq/raw/source/${BATCH}/fastqs/original
CLEANED_OUTPUT=/exports/igmm/eddie/Glioblastoma-WGS/RNA-seq/raw/source/${BATCH}/fastqs/cleaned
QC_DIR=/exports/igmm/eddie/Glioblastoma-WGS/RNA-seq/raw/source/${BATCH}/qc/fastp

echo "Detecting adaptor sequences and removing them..."
fastp -i ${READ_OUTPUT}/${PATIENT_ID}_R1.fastq.gz -I ${READ_OUTPUT}/${PATIENT_ID}_R2.fastq.gz \
     -o ${CLEANED_OUTPUT}/${PATIENT_ID}_R1.fastq.gz -O ${CLEANED_OUTPUT}/${PATIENT_ID}_R2.fastq.gz \
     -h ${QC_DIR}/${PATIENT_ID}.fastp.html

echo "Performing fastqc on cleaned fastqs..."
module load igmm/apps/FastQC/0.11.9
fastqc -t 20 ${CLEANED_OUTPUT}/${PATIENT_ID}_R1.fastq.gz
fastqc -t 20 ${CLEANED_OUTPUT}/${PATIENT_ID}_R2.fastq.gz
