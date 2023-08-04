#!/bin/bash

#This script uses various HOMER tools.
#To run this script, do qsub -t <1-n> submit_homer.sh ...

#$ -N homer_annotate
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
HOMER_DIR=$3
FASTA_DIR=$4

source $CONFIG

SAMPLE_ID=`head -n $SGE_TASK_ID $BED_LIST | cut -f 1 | tail -n 1`
BED_FILE=`head -n $SGE_TASK_ID $BED_LIST | cut -f 2 | tail -n 1`

#### Annotate (need to use unzipped fasta file for the reference, and unzipped gencode gtf file)
#### Based on http://homer.ucsd.edu/homer/ngs/annotation.html
annotatePeaks.pl $BED_FILE mm10 -CpG -annStats $HOMER_DIR/${SAMPLE_ID}.annotated.stats.txt  > $HOMER_DIR/${SAMPLE_ID}.annotated.txt

#### Convert to fasta for downstream use
bedtools getfasta -fi $MM10_REF_FASTA -bed $BED_FILE > $FASTA_DIR/${SAMPLE_ID}.fasta




