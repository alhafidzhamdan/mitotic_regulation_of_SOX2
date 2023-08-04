#!/bin/bash

# Attempt to run quantification of chip-seq reads within telomeric regions
# Uses telomerehunter based on https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2851-0#Sec9 and 
# https://www.nature.com/articles/s41467-019-13824-9#Sec11 (<- RNA seq, TERRA reads)

# To run this script, do 
# qsub -t 1-n submit_telomerehunter_chip_seq.sh IDS

#$ -N teloHunt_chip
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -l h_rt=64:00:00

export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/snakemake/bin:$PATH
export PATH=/exports/igmm/eddie/Glioblastoma-WGS/anaconda/envs/telomerehunter/bin:$PATH

IDS=$1

## Bam files are sorted unfiltered merged bams
SAMPLE_ID=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 1`
CHIP_BAM=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 2`
INPUT_BAM=`head -n $SGE_TASK_ID $IDS | tail -n 1 | cut -f 3`

## Based on article https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2851-0#Abs1
## And code/explanation is based here https://www.dkfz.de/en/applied-bioinformatics/telomerehunter/telomerehunter.html#section3
## TelomereHunter takes alignment information into account and reports the abundance of variant repeats in telomeric sequences. 
## Repeat threshold -rt set to 20 based on 140bp read length

OUTPUT_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/mitotic_bookmarking/telomerehunter
MM10_CYTOBAND=/exports/igmm/eddie/Glioblastoma-WGS/resources/refgenome_mmusculus/mm10_cytoBand.txt.gz

telomerehunter \
    -ibt $CHIP_BAM \
    -ibc $INPUT_BAM \
    -rt 20 \
    -o $OUTPUT_DIR \
    -p ${SAMPLE_ID} \
    -b $MM10_CYTOBAND




