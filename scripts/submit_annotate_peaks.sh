#!/bin/bash

#This script uses various HOMER tools.
#To run this script, do qsub -t <1-n> submit_annotate_peaks.sh ...

#$ -N HOMER_ANNOTATE
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -l h_vmem=4G
#$ -l h_rt=2:00:00
#$ -pe sharedmem 4

unset MODULEPATH
. /etc/profile.d/modules.sh

module load igmm/apps/homer/4.10

CONFIG=$1
BED_LIST=$2

source $CONFIG

SAMPLE_ID=`head -n $SGE_TASK_ID $BED_LIST | cut -f 1 | tail -n 1`
BED_FILE=`head -n $SGE_TASK_ID $BED_LIST | cut -f 2 | tail -n 1`
HOMER_DIR=/exports/igmm/eddie/Glioblastoma-WGS/ChIP-seq/SOX_TF/downstream/homer 

#### Annotate 
annotatePeaks.pl $BED_FILE hg38 -CpG -annStats $HOMER_DIR/${SAMPLE_ID}.annotated.stats.txt  > $HOMER_DIR/${SAMPLE_ID}.annotated.txt

#### Subset
cat $HOMER_DIR/${SAMPLE_ID}.annotated.txt | grep -si promoter | cut -f 2-4 | grep -vsi start > $HOMER_DIR/${SAMPLE_ID}.annotated.promoter.bed
cat $HOMER_DIR/${SAMPLE_ID}.annotated.txt | grep -si intergenic | cut -f 2-4 | grep -vsi start > $HOMER_DIR/${SAMPLE_ID}.annotated.intergenic.bed
cat $HOMER_DIR/${SAMPLE_ID}.annotated.txt | grep -vsi intergenic | grep -vsi promoter | cut -f 2-4 | grep -vsi start > $HOMER_DIR/${SAMPLE_ID}.annotated.intragenic.bed
    
#### Convert to fasta
bedtools getfasta -fi $REFERENCE -bed $HOMER_DIR/${SAMPLE_ID}.annotated.promoter.bed > $HOMER_DIR/${SAMPLE_ID}.annotated.promoter.fasta
bedtools getfasta -fi $REFERENCE -bed $HOMER_DIR/${SAMPLE_ID}.annotated.intergenic.bed > $HOMER_DIR/${SAMPLE_ID}.annotated.intergenic.fasta
bedtools getfasta -fi $REFERENCE -bed $HOMER_DIR/${SAMPLE_ID}.annotated.intragenic.bed > $HOMER_DIR/${SAMPLE_ID}.annotated.intragenic.fasta


